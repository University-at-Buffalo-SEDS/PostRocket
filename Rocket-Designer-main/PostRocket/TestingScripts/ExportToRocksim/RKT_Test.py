import xml.etree.ElementTree as ET



rocksim = ET.Element('RockSimDocument')
di = ET.SubElement(rocksim,'DesignInformation')
rd = ET.SubElement(di,'RocketDesign')
s3p = ET.SubElement(rd,'Stage3Parts')


