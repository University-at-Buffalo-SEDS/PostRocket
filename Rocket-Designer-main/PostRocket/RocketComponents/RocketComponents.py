# Authors: Alex Post, Jasper Stedma 7/15/2024
# The purpose of this file is to set up classes where a user can input a rocket design.

#######################################################################################################################
#######################################################################################################################
# Imports

import numpy as np
import scipy as sci
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.transforms as tfm
import pickle
import os
import xml.etree.ElementTree as ET
# from PostRocket.Propulsion.Funwithpropulsion import TNAM # Fake propulsion curve function for testing purposes
#from PostRocket.Propulsion.TNAM_v2func import TNAM # Proper TNAM function from the Magnani Modified McNalley Method (MMMM or M4)
# from PostRocket.Propulsion.HRAP_sim import TNAM # Python implementation of the HRAP matlab application
from PostRocket.Propulsion.HRAP_Matlab_sim import TNAM # Functional implementation of the HRAP matlab application
from PostRocket.FlightSim.Drag import *
from PostRocket.FlightSim.Aero import *

class material:
    def __init__(self,density,max_stress,surface_roughness):
        self.density = density
        self.max_stress = max_stress
        self.surface_roughness = surface_roughness


class sheet_material(material):
    def __init__(self,area_density,thickness,max_stress,surface_roughness):
        self.area_density = area_density
        density = area_density/thickness
        super(sheet_material,self).__init__(density,max_stress,surface_roughness)


#######################################################################################################################
#######################################################################################################################
# Definitions of classes for components inside the rocket, components of the airframe of the rocket, and the collective rocket
# Order of instantiation should be internal -> external -> rocket

class component:
    # Create a template for components
    def __init__(self,mass,cg,length,diameter):
        'mass [kg], cg [m from front], length [m], diameter [m]'
        self.mass = mass # Dry mass not including consumables like propellant [kg]
        self.cg = cg # cg relative to component's front [m]
        self.length = length # Component length [m]
        self.diameter = diameter # Component diameter (cylinder assumed for all) [m]
    
    def set_length(self,new_length):
        self.length = new_length
    
    def calc_total_mass(self):
        return self.mass
        
    
class internal_component(component): # Base internal component should be a simple mass object
    def __init__(self,mass,cg,length,diameter,*location):
        '''mass [kg], cg [m from front], length [m], diameter [m], location [m from parent front]

           Component that is inside the airframe. Most basic form is a miscellaneous point mass object'''
        super().__init__(mass,cg,length,diameter)
        self.Iyy = mass * ( (0.75*diameter**2 + length**2)/12) # Solid cylinder Iyy abt the top
        self.location = location[0] # Distance from the front of a parent part to the front of this part
        self.parent_component = [] # parent component is assigned when self is passed to an external component
    
    def install_in(self,parent,location):
        self.location = location
        parent.add_subcomponent(self)

    def get_draw_points(self, total_length):
        x = np.array([0, 0, self.length, self.length]) + total_length + self.location
        y = np.array([-self.diameter, self.diameter, self.diameter, -self.diameter]) / 2
        return x,y    


#### NOTE: Might want to redefine around CG placement instead of front

class external_component(component):
    def __init__(self,length,diameter,thickness,material,*subcomponents):
        '''length [m], diameter [m], thickness [m], material, subcomponent1, subcomponent2,...

           Component that is part of the airframe. Most basic form is essentially a body tube. Contains internal components.'''
        super().__init__([],[],length,diameter)
        self.thickness = thickness # Component wall thickness [mm]
        self.material = material
        self.calc_mass()
        self.cg = length/2
        self.subcomponents = [] # Initializes subcomponent list to append all passed
        self.add_subcomponents(*subcomponents)
        self.calc_total_mass()
        self.calc_total_cg()
    
    def calc_mass(self):
        self.mass = np.pi * self.diameter * self.length * self.thickness * self.material.density
        
    def calc_total_mass(self):
        total_mass = self.mass
        for subcomponent in self.subcomponents:
            total_mass += subcomponent.mass
        self.total_mass = total_mass
        return total_mass

    def calc_total_cg(self):
        total_mcg = self.mass*self.cg
        for subcomponent in self.subcomponents:
            total_mcg += subcomponent.mass*(subcomponent.cg + subcomponent.location)
        self.total_cg = total_mcg / self.total_mass
        return self.total_cg
    
    def add_subcomponents(self,*subcomponents):
        for subcomponent in subcomponents:
            if subcomponent.parent_component:
                raise Exception('Subcomponent already assigned to parent component!')
            else:
                self.subcomponents.append(subcomponent)
                subcomponent.parent_component = self # Assigns self as the parent of each subcomponent
    
    def get_draw_points(self, total_length):
        x = np.array([0, 0, self.length, self.length]) + total_length
        y = np.array([-self.diameter, self.diameter, self.diameter, -self.diameter])/2
        return x,y

    
    # CP function might need AoA input depending on whether those relations are computed
    # for the rocket or for individual components  


#######################################################################################################################
#######################################################################################################################
# Definitions of aerodynamic components

class von_karman():
    '''Contains functions to calculate properties of a Von-Karman (LD Haack)'''
    def radius(R,L,x):
        '''Von Karman curve equation for radius at distance x from the nose give base radius R and total length L'''
        theta = np.arccos( 1 - 2*x/L )
        Y = R/np.sqrt(np.pi) * np.sqrt( theta - 0.5*np.sin(2*theta) ) # radius
        return Y
    
    def slope(R,L,x):
        '''Derivative of the Von Karman curve at distance x from the nose give base radius R and total length L'''
        theta = np.arccos( 1 - 2*x/L )
        dydx = -R*( np.cos(2*theta) - 1 ) / ( np.sqrt( (2*np.pi) * (x*(L-x)/(L**2)) * (2*theta - np.sin(2*theta)) )*L )
        return dydx


    def surface_area(R,L_tip,L):
        '''Surface area of a (potentially) truncated Von Karman curve with a base radius R, total length L, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/1000 # Step for x

        X = np.arange(L_tip,L,dx) # x coordinates
        Y = np.zeros(len(X)) # nosecone radius initialized
        dYdX = np.zeros(len(X)) # radius derivatives initialized
        S = 0 # integral surface area

        for x in X:
            Y = von_karman.radius(R,L,x) # radius 
            dYdX = von_karman.slope(R,L,x) # Derivative
            dS = 2*np.pi*Y*np.sqrt( 1 + dYdX**2 ) * dx
            S += dS
        return S
    
    def plane_area(R,L_tip,L):
        '''Planiform area of a (potentially) truncated Von Karman curve with a base radius R, total length L, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/1000 # Step for x

        X = np.arange(0,L,dx) # x coordinates
        A = 0 # integral area

        for x in X:
            Y = von_karman.radius(R,L,x)
            dA = Y * dx
            A += dA
        return A
    
    def volume(R,L_tip,L):
        '''Total volume enclosed by a (potentially) truncated Von Karman curve with a base radius R, total length L, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/1000 # Step for x

        X = np.arange(L_tip,L,dx) # x coordinates
        V = 0 # integral volume

        for x in X:
            Y = von_karman.radius(R,L,x)
            dV = np.pi * Y**2 * dx
            V += dV
        return V
    
    def Iyy_from_tip(R,L_tip,L,t):
        '''Moment of inertia / density of a (potentially) truncated Von Karman curve relative to its total tip length with a base radius R, total length L, wall thickness t, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/1000 # Step for x

        X = np.arange(L_tip,L,dx) # x coordinates
        I = 0 # integral volume
        for x in X:
            Y = von_karman.radius(R,L,x)
            if t == 0:
                dV = np.pi * Y**2 * dx
            else:
                dV = np.pi * (2*Y*t - t**2 ) * dx

            I += dV*x*x

        return I
    

    
    def cg(R,L_tip,L,t):
        '''Center of gravity and moment of inertia of a (potentially) truncated Von Karman curve with a base radius R, total length L, wall thickness t, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/1000 # Step for x
        X = np.arange(L_tip,L,dx) # x coordinates
        V = 0 # integral volume
        VX = 0 # integral volume times distance from tip
        for x in X:
            Y = von_karman.radius(R,L,x)
            if t == 0:
                dV = np.pi * Y**2 * dx
            else:
                dV = np.pi * (2*Y*t - t**2 ) * dx

            V += dV
            VX += dV*x

        return VX / V
    
    
    def cp(R,L_tip,L):
        '''Center of planiform area approximation of the center of pressure for a (potentially) truncated Von Karman curve with a base radius R, total length L, and length of the removed portion L_tip from 0 to L'''
        if 0 > L_tip or L_tip > L:
            return -1
        
        dx = L/10 # Step for x
        X = np.arange(L_tip,L,dx) # x coordinates
        A = 0 # integral area
        AX = 0 # integral area times distance from tip
        for x in X:
            Y = von_karman.radius(R,L,x)
            dA = Y * dx
            A += dA
            AX += dA*x

        return AX / A


class tipped_nosecone(external_component):
    def __init__(self,length,diameter,thickness,material,tip_length,tip_material,*subcomponents):
        '''length [m], base diameter [m], wall thickness [m], wall material, tip length [m], tip material, subcomponent1,... subcomponentN
        
            Hollow Von Karman nosecone with a solid tip of a different material, usually metal'''
        super().__init__(length,diameter,thickness,material,*subcomponents)
        # Base shell values
        base_mass = material.area_density * von_karman.surface_area(diameter/2,tip_length,length)
        base_cg = von_karman.cg(diameter/2,tip_length,length,thickness) + tip_length
        base_Iyy = material.density * von_karman.Iyy_from_tip(diameter/2,tip_length,length,thickness) - base_mass*base_cg**2
        # Nosecone tip values
        tip_radius = von_karman.radius(diameter/2,length,tip_length)
        tip_mass = tip_material.density * von_karman.volume(tip_radius,0,tip_length)
        tip_cg = von_karman.cg(tip_radius,0,tip_length,0)
        tip_Iyy = tip_material.density * von_karman.Iyy_from_tip(tip_radius,0,tip_length,0) - tip_mass*tip_cg**2
        # Total values
        self.mass = tip_mass + base_mass
        self.cg = (tip_mass*tip_cg + base_mass*base_cg) / self.mass
        self.Iyy = (base_Iyy + base_mass*(self.cg - base_cg)**2) + (tip_Iyy + tip_mass*(self.cg - tip_cg)**2)
        # Aerodynamic properties
        self.reference_area = np.pi * diameter**2 / 4
        self.cp = von_karman.cp(diameter/2,0,length)
        plane_area = von_karman.plane_area(diameter/2,0,length)
        self.cA = self.cp # center of plane area
        self.ArefCNa, self.ArefCNa2 = CNderivatives_nosecone(self.reference_area, plane_area, diameter)
        

    def get_draw_points(self, total_length):
        numPts = 50
        x = np.linspace(0,self.length,numPts)
        y = von_karman.radius(self.diameter/2, self.length, x)
        x = np.append( x, np.flip(x) )
        y = np.append( y,-np.flip(y) )
        x += total_length
        return x,y

        
class body_tube(external_component):
    def __init__(self,length,diameter,thickness,material,*subcomponents:component):
        super(body_tube,self).__init__(length,diameter,thickness,material,*subcomponents)
        self.set_length(self.length)

    def calc_Iyy(self):
        self.Iyy = self.mass*( 3*(0.5*self.diameter**2 - self.diameter*self.thickness + self.thickness**2) + self.length**2 )/12

    def set_length(self,new_length):
        self.cg = new_length/2
        self.length = new_length
        self.plane_area = self.length * self.diameter
        self.calc_mass()
        self.calc_Iyy()
        # Aerodynamic properties
        self.cp = self.length / 2
        self.cA = self.cp # center of plane area
        self.ArefCNa, self.ArefCNa2 = CNderivatives_bodytube(self.plane_area)     


class von_karman_boattail(external_component):
    def __init__(self,length,fore_diameter,aft_diameter,thickness,material,*subcomponents:component):
        '''length [m], fore diameter [m], aft diameter [m], wall thickness [m], wall material, subcomponent1, subcomponent2,...

           Von Karman boattail with the same parameters as a similarly-sized conical boattail'''
        super().__init__(length,fore_diameter,thickness,material,*subcomponents)
        self.fore_diameter = fore_diameter
        self.aft_diameter = aft_diameter
        von_karman_length = self.find_boattail_vk_length(length,fore_diameter,aft_diameter)
        self.von_karman_length = von_karman_length


        self.mass = material.density * thickness * von_karman.surface_area(fore_diameter/2,von_karman_length - length, von_karman_length)
        cg_from_tip = von_karman.cg(fore_diameter/2,von_karman_length - length,von_karman_length,thickness) # cg of the boattail measured from the theoretical tip of its full von karman curve
        self.cg =  von_karman_length - cg_from_tip # Actual cg from the boattail front 
        Iyy_from_tip = material.density * von_karman.Iyy_from_tip(fore_diameter/2,von_karman_length - length,von_karman_length,thickness) # Iyy measured from the theoretical tip of its full von karman curve
        self.Iyy = Iyy_from_tip - self.mass*cg_from_tip**2 # Iyy from the proper cg
        # Aerodynamic properties
        diameter_ratio = fore_diameter / aft_diameter
        self.cp = length/3 * ( 1 + ( 1 - diameter_ratio )/( 1 - diameter_ratio**2 ) )
        plane_area = von_karman.plane_area(fore_diameter/2,von_karman_length - length, von_karman_length)
        self.cA = von_karman.cp(fore_diameter/2,von_karman_length - length, von_karman_length) # center of plane area
        self.ArefCNa, self.ArefCNa2 = CNderivatives_boattail(plane_area, diameter_ratio, fore_diameter)
        

    def find_boattail_vk_length(self,L_bt,d1,d2):
        '''Finds the total length of a von karman curve with base diameter d1 that, when truncated to desired length L_bt, has an aft diameter d2'''
        r1 = d1/2
        r2 = d2/2
        L_vk = L_bt # von karman length required to get desired boattail dimensions
        r_vk = 0 # radius at the desired boattail length
        dL_vk = L_bt/1000 # step size for the total von karman length

        while r_vk < r2:
            L_vk += dL_vk
            r_vk = von_karman.radius(r1,L_vk,L_vk-L_bt)
        
        return L_vk
    
    def get_draw_points(self, total_length):
        numPts = 50
        x = np.linspace(0,self.length,numPts)
        y = von_karman.radius(self.fore_diameter/2, self.von_karman_length, self.von_karman_length - x )
        x = np.append( x, np.flip(x) )
        y = np.append( y,-np.flip(y) )
        x += total_length
        return x,y

    
class composite_fin_group(external_component):
    def __init__(self,number,height,root_chord,tip_chord,sweep_length,thickness,skin_material,core_material,parent_diameter,distance_from_base):
        '''Number of fins, height/semispan [m], root chord [m], tip chord [m], sweep length [m from root to tip leading edge], thickness [m], material
           Trapezoidal fin set with basic dimensions'''
        self.number = number
        self.thickness = thickness
        self.height = height
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.sweep_length = sweep_length
        self.skin_material = skin_material
        self.core_material = core_material
        self.parent_radius = parent_diameter/2
        self.distance_from_base = distance_from_base
        self.calc_midchord()
        self.calc_planform_area()
        self.calc_cp(self.height,self.root_chord,self.tip_chord,self.sweep_length,self.parent_radius)
        self.calc_Iyy()
        self.calc_Aero()
        
        super().__init__(0,0,thickness,skin_material) # zero length in order to be placed on the previous component when passed to a rocket

    def calc_total_mass(self):
        self.mass = (self.root_chord + self.tip_chord)*self.height/2 * self.thickness * self.material.density
        return self.mass
    
    def calc_cp(self,s,a,b,m,R): # all barrowman for now
        self.cA = m/3 * (a + 2*b)/(a+b) + 1/6 * (a + b - a*b/(a+b))
        self.cp = self.cA - a - self.distance_from_base # subtracts root chord and base distance at the end to account for being at the bottom of the body tube 
    
    def calc_total_cg(self):
        self.cg = self.cA - self.root_chord - self.distance_from_base
        return self.cg

    def calc_Iyy(self):
        x = np.array([0, self.sweep_length, self.sweep_length + self.tip_chord, self.root_chord])
        y = np.array([0, self.height, self.height, 0])
        dx = (max(x) - min(x))/1000
        Iyy = 0 # integral area Iyy (rocket's y axis align with the fin's z axis for this calculation)
        for x_val in np.arange(min(x),max(x)-dx,dx):
            y_val0 = np.interp(x_val,x,y)
            y_val1 = np.interp(x_val + dx,x,y)
            dA = 0.5*(y_val0 + y_val1) * dx
            Iyy += dA * (x_val + dx/2 - self.cA)**2
        
        self.Iyy = Iyy * self.thickness * self.skin_material.density * self.number

    def calc_Aero(self):
        self.ArefCNa, self.ArefCNa2 = CNderivatives_fins(self.planform_area, self.number, self.height, self.sweep_length, self.root_chord, self.tip_chord, self.parent_radius)
        

    def get_drag(self, velocity, density, temperature, Rs, rocket_length):
        cd_fin, cd_i = fin_drag(self.thickness, (self.root_chord + self.tip_chord)/2, self.number, self.planform_area, self.exposed_area, self.parent_radius*2, Rs, velocity, density, temperature, rocket_length)
        return cd_fin, cd_i

    def stability_error(self, semispan, Cnareq):
        ArefCNa, ArefCNa2 = CNderivatives_fins(self.planform_area, self.number, semispan, self.sweep_length, self.root_chord, self.tip_chord, self.parent_radius)
        error = ArefCNa - Cnareq
        return error

    def update_height(self, Cnareq):
        matched_height = root_scalar(self.stability_error,args=(Cnareq), x0=self.height)
        self.height = matched_height.root
        self.calc_Aero()
        self.calc_midchord()
    
    def calc_midchord(self):
        self.midchord = np.sqrt((self.sweep_length + (self.tip_chord - self.root_chord)/2)**2 + self.height**2)

    def calc_planform_area(self):
        self.planform_area = 0.5 * (self.root_chord + self.tip_chord) * self.height

    # TODO: Determine which dimension will be optimized to maintain Cp (likely fin height)
    # TODO: Implement Jaynie's code to compute mass from the Nomex and Fin_CF materials
    # TODO: Possibly ACTUALLY implement mass scaling based on area scaling

    def get_draw_points(self, total_length):
        
        x = np.array([0, self.sweep_length, self.sweep_length + self.tip_chord, self.root_chord]) + total_length - self.root_chord - self.distance_from_base
        y = np.array([self.parent_radius, self.parent_radius + self.height, self.parent_radius + self.height, self.parent_radius])

        x0 = np.array(x)
        y0 = np.array(y)

        num_fins_visible = self.number // 2
        angle_between_fins = np.pi * 2 / self.number
        for i in range(num_fins_visible):
            x = np.append( x, x0)
            y = np.append( y, y0 * np.cos(angle_between_fins * (i+1)) )
            
        return x,y


class parachute(internal_component):
    def __init__(self,packed_length,packed_diameter,location,material,cd):

        self.cd = cd
        self.material = material
        self.landing_speed = 7.5 # [m/s]
        mass = 0 # Starting mass of zero because total rocket mass is unknown
        super().__init__(mass,packed_length/2,packed_length,packed_diameter,location)
        
    
    def calc_total_mass(self,rocket_mass,landing_speed):
        g = 9.81 # gravity [m/s^2]
        rho = 1.293 # rough air density [kg/m^3]
        area_1 = 2*g/(self.cd*rho*landing_speed**2)
        area = area_1 / (1 - area_1*self.material.area_density) * rocket_mass # Area calculation accounts for the parachute's own increased mass with the 1/(1-a*rho) term
        mass = area * self.material.area_density
        self.mass = mass
        return mass
    

    
#######################################################################################################################
#######################################################################################################################
# Definitions of propulsion classes

class propulsion_element():
    def __init__(self,length,diameter,material):
        ''' Propulsion elements take an evenly spaced massurve vector and burn time and calculate their current mass based on the time elapsed since ignition'''
        self.length = length
        self.diameter = diameter
        self.material = material
        self.pressure = []
        self.mass_curve = []
        self.burn_time = []
        self.dt = []

    def calc_thickness(self):
        self.thickness = 2 * self.pressure * self.diameter / (2 * self.material.max_stress) # hoop stress determines thickness with a 2x factor of safety
    
    def calc_mass(self):
        mass =  np.pi*(self.diameter*self.length + (self.diameter**2)/2 )*self.thickness*self.material.density # Contained, thin-walled cylinder
        self.mass = mass

    def calc_internal_area(): # different for every type of motor and geometry
        return
        
    def update_element(self,mass_curve,pressure_curve,burn_time):
        self.mass_curve = mass_curve
        self.burn_time = burn_time        
        self.pressure_curve = pressure_curve
        self.pressure = max(pressure_curve)

    def interp_curve(self,t,curve): 
        value = np.interp(t,self.burn_time,curve)
        return value
        
    def get_propellant_mass(self,t):
        mass = self.interp_curve(t,self.mass_curve)
        return mass

    def calc_Iyy(self):
        wall_Iyy = self.mass*( 3*(0.5*self.diameter**2 - self.diameter*self.thickness + self.thickness**2) + self.length**2 )/12
        bulkheads_Iyy = 2*self.mass*( ( 0.75*self.diameter**2 + self.thickness**2 )/12 + (self.length**2)/4 )
        self.Iyy = wall_Iyy + bulkheads_Iyy

class oxidizer_tank(propulsion_element,body_tube):
    def __init__(self,length,diameter,material):
        '''Length [m], outer diameter [m], wall material, oxidizer propellant'''
        super().__init__(length,diameter,material)
        self.subcomponents = []
        
    
    def calc_internal_area(self):
        area = np.pi * (self.diameter/2 - self.thickness)**2
        return area

    def update_element(self,mass_curve,pressure_curve,burn_time,oxidizer_volume):
        super().update_element(mass_curve,pressure_curve,burn_time)
        
        # Setting dimensions
        self.calc_thickness()
        self.length = oxidizer_volume/ self.calc_internal_area()
        self.set_length(self.length)
        self.calc_mass()
        self.calc_Iyy()
        self.cg = self.length / 2


class combustion_chamber(propulsion_element,internal_component):
    def __init__(self,length,diameter,prechamber_length,postchamber_length,material):
        '''Length [m], grain diameter [m], pre-chamber length [m], post-chamber length [m], wall material, fuel propellant'''
        cg = length/2
        internal_component.__init__(self,0,cg,length,diameter,0) # location 0 because chamber extends through its whole body tube, mass 0 until initialized
        propulsion_element.__init__(self,length,diameter,material)
        self.thickness = []
        self.grain_length = []
        self.port_diameter = []
        self.grain_cg = []
        self.prechamber_length = prechamber_length
        self.postchamber_length = postchamber_length

    def update_element(self,mass_curve,burn_time,pressure_curve,grain_length,thrust_curve,port_diameter):
        super().update_element(mass_curve,pressure_curve,burn_time)
        self.grain_length = grain_length
        self.port_diameter = port_diameter
        self.length = grain_length + self.prechamber_length + self.postchamber_length
        self.thrust_curve = thrust_curve
        self.calc_thickness()
        self.calc_mass()
        self.calc_Iyy()
        self.cg = self.length/2
        self.grain_cg = self.prechamber_length + grain_length/2
        self.parent_component.set_length(self.length)

    def get_pressure(self,t):
        pressure = self.interp_curve(t,self.pressure_curve)
        return pressure
    
    def get_thrust(self,t):
        thrust = self.interp_curve(t,self.thrust_curve)
        return thrust


class nozzle(internal_component):
    def __init__(self,throat_area,exit_area,mass,length,diameter,*location):
        super().__init__(mass,length/2,length,diameter,*location)
        self.throat_area = throat_area
        self.exit_area = exit_area
        
        
#######################################################################################################################
#######################################################################################################################
# Definition of the rocket class

class rocket:
    def __init__(self, *components):
        '''component1, component2, ...

            Rocket is a container for all components and has methods for analyzing the whole structure''' 
        
        self.name = 'TNAM'
        self.components = [] # Initializes subcomponent list and appends all passed
        self.diameter = components[0].diameter
        self.reference_area = components[0].reference_area
        for component in components: # Iterating through every component from nose to tail
            self.components.append(component)
            if type(component) == oxidizer_tank:
                self.oxidizer_tank = component
            elif type(component) == tipped_nosecone:
                self.nosecone = component
            elif type(component) == von_karman_boattail:
                self.boattail = component
            elif type(component) == composite_fin_group:
                self.fins = component
            for subcomponent in component.subcomponents:
                if type(subcomponent) == combustion_chamber:
                        self.combustion_chamber = subcomponent
                if type(subcomponent) == parachute:
                        self.parachute = subcomponent
                if type(subcomponent) == nozzle:
                        self.nozzle = subcomponent

        self.mass = []
        self.cg = []
        self.oxidizer_cg = []
        self.fuel_cg = []
        self.mass = []
        self.length = []
        self.trajectory = []
        self.updated = False
        
    def calc_total_mass(self):
        total_mass = 0
        for component in self.components:
            if component is not self.parachute: # Adds total mass not including the parachute first
                total_mass += component.calc_total_mass()
        self.mass = total_mass + self.parachute.calc_total_mass(total_mass,6) # At the end, recalculates the parachute to land at 6 m/s   

    def calc_cg(self):
        total_mcg = 0
        total_length = 0
        for component in self.components:
            if component is self.oxidizer_tank:
                self.oxidizer_cg = total_length + component.cg
            else:
                for subcomponent in component.subcomponents:
                    if subcomponent is self.combustion_chamber:
                        self.fuel_cg = total_length + subcomponent.location + subcomponent.grain_cg

            total_mcg += component.calc_total_mass()*(component.calc_total_cg() + total_length)
            total_length += component.length
        self.cg = total_mcg / self.mass
    
    def calc_Iyy(self):
        total_Iyy = 0
        total_length = 0        
        for component in self.components:
            total_Iyy += component.Iyy + component.mass*(total_length + component.cg - self.cg)**2 # Adding component Iyy abt its cg 

            for subcomponent in component.subcomponents:

                total_Iyy += subcomponent.Iyy + subcomponent.mass*(total_length + subcomponent.location + subcomponent.cg - self.cg)**2 # adding subcomponent's Iyy abt its cg

            total_length += component.length
        
        self.Iyy = total_Iyy
    
    def calc_total_length(self):
        total_length = 0
        for component in self.components:
            total_length += component.length
        self.length = total_length
    
    def get_aero(self,t,q,AoA, mach, velocity, density, temperature, Rs):
        
        # Compressibility with capped prandtl correction factor 
        Prandtl_correction = min(5,abs(1-mach**2)**-0.5)
        # if (0.8 < mach) and (mach < 1.1):
        #     Prandtl_correction = (1-0.8**2)**-0.5
        # else:
        #     Prandtl_correction = abs(1-mach**2)**-0.5
        
        Drag = cd(self, t, q, mach,  velocity, density, temperature, Rs, AoA) * self.reference_area * q * Prandtl_correction

        if abs(AoA) < 1e-12: # If AoA goes to zero, default to Barrowman static Cp where AoA is factored out and body lift is not considered
            AoA = 1
            AoA_sq = 0
        else:
            AoA_sq = AoA * abs(AoA)

        # Lift and CP modeling
        total_length = 0
        total_ArefCN = 0
        total_ArefCN_Xcp = 0
        for component in self.components:
            ArefCN = (component.ArefCNa * AoA + component.ArefCNa2 * AoA_sq)
            ArefCN_Xcp = ArefCN * (total_length + component.cp)
            total_ArefCN += ArefCN      
            total_ArefCN_Xcp += ArefCN_Xcp
            total_length += component.length
        Lift = q * total_ArefCN * Prandtl_correction
        Cp = total_ArefCN_Xcp / total_ArefCN

        # Correcting Drag with the calculated lift
        Drag = (Drag - Lift*np.sin(AoA))/np.cos(AoA)

        return Lift, Drag, Cp
    
    def Fin_CNaMatch(self, stab):
        total_length = 0
        total_ArefCN = 0
        total_ArefCN_Xcp = 0
        for component in self.components:
            if component == self.components[-2]:
                ArefCNa_fin = component.ArefCNa
                Xcpfin = total_length + component.cp
                ArefCN_Xcp_fins = ArefCNa_fin * (total_length + component.cp)
            ArefCN = component.ArefCNa
            ArefCN_Xcp = ArefCN * (total_length + component.cp)
            total_ArefCN += ArefCN      
            total_ArefCN_Xcp += ArefCN_Xcp
            total_length += component.length
        left_over_ArefCN_Xcp = total_ArefCN_Xcp - ArefCN_Xcp_fins
        left_over_ArefCN = total_ArefCN - ArefCNa_fin
        self.initial_cp = total_ArefCN_Xcp / total_ArefCN # TEMP
        cg = self.get_cg(0)
        desired_cp = cg + stab*self.diameter
        self.initial_stability = (self.initial_cp - cg) / self.diameter # TEMP
        CN_desired = (desired_cp * left_over_ArefCN - left_over_ArefCN_Xcp ) / (Xcpfin - desired_cp)
        self.components[-2].update_height(CN_desired)
    
    def get_thrust(self,t,pressure):
        base_thrust = self.combustion_chamber.get_thrust(t)
        thrust_adjustment = self.nozzle.exit_area * (self.combustion_chamber.get_pressure(t) - pressure)
        total_thrust = base_thrust + thrust_adjustment
        return total_thrust
    
    def get_mass(self,t):
        fuel_mass = self.combustion_chamber.get_propellant_mass(t)
        oxidizer_mass = self.oxidizer_tank.get_propellant_mass(t)
        
        return self.mass + fuel_mass + oxidizer_mass

    def get_cg(self,t):
        fuel_mass = self.combustion_chamber.get_propellant_mass(t)
        oxidizer_mass = self.oxidizer_tank.get_propellant_mass(t)

        cg = (self.mass*self.cg + fuel_mass*self.fuel_cg + oxidizer_mass*self.oxidizer_cg) / (self.mass + fuel_mass + oxidizer_mass)

        return cg

    def get_Iyy(self,t):

         #self.get_mass(t) * (0.75*self.diameter**2 + self.length**2)/12  # Cylindrical approximation

        fuel_mass = self.combustion_chamber.get_propellant_mass(t)
        oxidizer_mass = self.oxidizer_tank.get_propellant_mass(t)
        current_cg = self.get_cg(t)

        # Oxidizer Contribution + cg adjustment for tank emptying
        oxidizer_radius = self.oxidizer_tank.diameter/2 - self.oxidizer_tank.thickness
        oxidizer_height = self.oxidizer_tank.length * oxidizer_mass / self.oxidizer_tank.get_propellant_mass(0)
        oxidizer_Iyy = oxidizer_mass * ( ( 3*oxidizer_radius**2 + oxidizer_height**2)/12 + (self.oxidizer_cg + (self.oxidizer_tank.length - oxidizer_height)/2 - current_cg)**2 )

        # Fuel Contribution + port radius adjustment for fuel burning
        fuel_outer_radius = self.combustion_chamber.diameter/2
        fuel_outer_area = np.pi*fuel_outer_radius**2
        initial_port_area = 0.25*np.pi*self.combustion_chamber.port_diameter**2
        fuel_area_remaining = (fuel_outer_area - initial_port_area) * fuel_mass / self.combustion_chamber.get_propellant_mass(0)
        fuel_inner_radius = np.sqrt( (fuel_outer_area - fuel_area_remaining)/np.pi )
        fuel_height = self.combustion_chamber.length
        fuel_Iyy = fuel_mass * ( ( 3*(fuel_inner_radius**2 + fuel_outer_radius**2) + fuel_height**2)/12 + (self.fuel_cg - current_cg)**2 )

        return self.Iyy + oxidizer_Iyy + fuel_Iyy

    def update_components(self, port_diameter,grain_length,oxidizer_volume):
        self.updated = True
        burn_time, oxidizer_pressure_curve, chamber_pressure_curve, oxidizer_mass_curve, fuel_mass_curve, thrust_curve = TNAM(  self.combustion_chamber.diameter , port_diameter , grain_length , oxidizer_volume  )
            ## TODO: ADD NOZZLE INPUTS TO TNAM PLS
        self.burn_time = burn_time[-1]
        self.burn_time_curve = burn_time
        self.oxidizer_tank.update_element(oxidizer_mass_curve,oxidizer_pressure_curve,burn_time,oxidizer_volume)
        self.combustion_chamber.update_element(fuel_mass_curve,burn_time,chamber_pressure_curve,grain_length,thrust_curve,port_diameter)
        self.calc_total_mass()
        self.calc_cg()
        self.calc_Iyy()
        self.calc_total_length()
        stability = 2 # Desired static stabiltiy at launch
        self.Fin_CNaMatch(stability)
        self.name = 'TNAM_' + str(round(port_diameter*1000,1)) + 'mm_'  + str(round(grain_length*1000)) + 'mm_'  + str(np.round(oxidizer_volume*1000)) + 'L'
    
    def draw(self, ax = 0):
        if ax == 0:
            ax = plt.axes()
        ax.axis('equal')
        ax.grid(True,'major',alpha=0.25)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True,'minor', alpha= 0.125)

        total_length = 0
        for component in self.components:
            x,y = component.get_draw_points(total_length)
            if type(component) == composite_fin_group: # Fin groups need to be split up to not have intersecting lines between fins
                for i in np.array(range(len(x)//4)) * 4:
                    xi = np.append(x[i:i+4],x[i])
                    yi = np.append(y[i:i+4],y[i])
                    ax.plot(xi, yi, linewidth=2, color='black')
            else:
                x = np.append(x,x[0])
                y = np.append(y,y[0])
                ax.plot(x, y, linewidth=2, color='black')

            if component is self.oxidizer_tank:
                ax.fill(x, y,color='blue',alpha=0.5)

            for subcomponent in component.subcomponents:
                x,y = subcomponent.get_draw_points(total_length)
                x = np.append(x,x[0])
                y = np.append(y,y[0])
                ax.plot(x, y, linewidth=1, color='black', alpha = 0.5)
                if subcomponent is self.combustion_chamber:
                    x[[0,1,4]] += subcomponent.prechamber_length
                    x[[2,3]] -= subcomponent.postchamber_length        
                    ax.fill(x, y,color='orange',alpha=0.5)
                    y[[0,3,4]] = -subcomponent.port_diameter / 2
                    y[[1,2]] = subcomponent.port_diameter / 2      
                    ax.fill(x, y,color='orange',alpha=0.55)
                
                if subcomponent is self.parachute:
                    ax.fill(x, y,color='black',alpha=0.25)
            
            total_length += component.length
            ax.text(self.length/2, self.diameter*2, self.name, horizontalalignment='center', fontsize=20, color='green')
        return ax

    def animate(self, ax = plt.axes()):
       
            ax.axis('equal')
            ax.grid(True,'major',alpha=0.25)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.grid(True,'minor', alpha= 0.125)
            fills = []
            total_length = 0
            for component in self.components:
                x,y = component.get_draw_points(total_length)
                if type(component) == composite_fin_group: # Fin groups need to be split up to not have intersecting lines between fins
                    for i in np.array(range(len(x)//4)) * 4:
                        xi = np.append(x[i:i+4],x[i])
                        yi = np.append(y[i:i+4],y[i])
                        ax.plot(xi, yi, linewidth=2, color='black')
                else:
                    x = np.append(x,x[0])
                    y = np.append(y,y[0])
                    ax.plot(x, y, linewidth=2, color='black')

                if component is self.oxidizer_tank:
                    fills.append(ax.fill(x, y,color='blue',alpha=0.5))

                for subcomponent in component.subcomponents:
                    x,y = subcomponent.get_draw_points(total_length)
                    x = np.append(x,x[0])
                    y = np.append(y,y[0])
                    ax.plot(x, y, linewidth=1, color='black', alpha = 0.5)
                    if subcomponent is self.combustion_chamber:
                        x[[0,1,4]] += subcomponent.prechamber_length
                        x[[2,3]] -= subcomponent.postchamber_length        
                        ax.fill(x, y,color='orange',alpha=0.5)
                        y[[0,3,4]] = -subcomponent.port_diameter / 2
                        y[[1,2]] = subcomponent.port_diameter / 2      
                        fills.append(ax.fill(x, y,color='orange',alpha=0.55))
                    
                    if subcomponent is self.parachute:
                        fills.append(ax.fill(x, y,color='black',alpha=0.25))
                
                total_length += component.length

                
            return ax
    
    def save(self, name = None, directory = None, trajectory = []):
        '''Save Rocket object as a serialized .obj file'''
        if name == None:
            name = self.name
        if directory == None:
            directory = "PostRocket/SavedRockets"

        export_folder  = f"{directory}/"

        if not os.path.exists(export_folder):
            os.makedirs(export_folder)
        
        # pickling the Rocket object
        self.trajectory = trajectory
        pickle.dump(self, open(f"{export_folder}/{name}.obj","wb") )

        if not self.updated:
            print("Warning: Rocket un-initialized, no motor data saved!")
    
    def load(directory):
        return pickle.load(open(directory,"rb"))

    def export(self, name = None, directory = None, trajectory = None):
        '''Export Rocket as a RockSim-compatible .rkt and motor as .rse'''
        if name == None:
            name = self.name
        if directory == None:
            directory = "PostRocket/SavedRockets"
        
        export_folder  = f"{directory}/{name}/"

        if not os.path.exists(export_folder):
            os.makedirs(export_folder)
            ## TODO: ADD ALL THE XML STUFF TMM
        
        rocksim = ET.Element('RockSimDocument')
        di = ET.SubElement(rocksim,'DesignInformation')
        rd = ET.SubElement(di,'RocketDesign')
        s3p = ET.SubElement(rd,'Stage3Parts')
        
        def newSE(element,field,value):
            ET.SubElement(element,field).text = str(round(value,3))
        
        for i,component in enumerate(self.components):
            if type(component) == tipped_nosecone:
                nose = ET.SubElement(s3p,'NoseCone')
                newSE(nose,"Len",component.length)
                newSE(nose,"BaseDia",component.diameter)
                newSE(nose,"ShapeCode",6)
            elif type(component) == body_tube:
                tube = ET.SubElement(s3p,"BodyTube")
                newSE(tube,"Len",component.length)
                newSE(tube,"OD",component.diameter)
                if self.combustion_chamber in component.subcomponents:
                    newSE(tube,"IsMotorMount",1)
            elif type(component) == oxidizer_tank:
                ox_tank = ET.SubElement(s3p,"BodyTube")
                newSE(ox_tank,"Len",component.length)
                newSE(ox_tank,"OD",component.diameter)
                newSE(ox_tank,"IsMotorMount",1)
            elif type(component) == composite_fin_group:
                ap = ET.SubElement(tube,"AttachedParts")
                fins = ET.SubElement(ap,"FinSet")
                newSE(fins,"FinCount",component.number)
                newSE(fins,"RootChord",component.root_chord)
                newSE(fins,"TipChord",component.tip_chord)
                newSE(fins,"SemiSpan",component.height)
                newSE(fins,"SweepDistance",component.sweep_length)
                newSE(fins,"Thickness",component.thickness)
                newSE(fins,"Xb",component.distance_from_base)
                newSE(fins,"LocationMode",2)
            elif type(component) == von_karman_boattail:
                tail = ET.SubElement(s3p,'NoseCone')
                newSE(tail,"Len",component.length)
                newSE(tail,"FrontDia",component.fore_diameter)
                newSE(tail,"RearDia",component.aft_diameter)
                newSE(tail,"ShapeCode",6)
            








            

