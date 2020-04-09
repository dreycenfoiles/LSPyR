from netgen.csg import *
from ngsolve import Mesh
from Materials import *


#Rename nanoparticle variables if more materials are implemented 
#Could maybe be improved by changing water and PMl domains into cylinders to save on space

def Nanosphere(particle):

	physical_space = particle + 150
	domain = physical_space + 100

	geo = CSGeometry()

	sphere1 = Sphere(Pnt(0,0,0),particle)
	sphere2 = Sphere(Pnt(0,0,0),physical_space)
	sphere3 = Sphere(Pnt(0,0,0),domain).bc('outer')

	AuNP = sphere1.mat('gold')
	water = (sphere2 - sphere1).mat('water')
	pml = (sphere3 - sphere2).mat('pml').maxh(80)

	geo.Add(AuNP)
	geo.Add(water)
	geo.Add(pml)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)

	return mesh


def Nanorod(aeff,ratio,mt_length):

	geo = CSGeometry()

	if mt_length == 0:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		physical_space = length + 100
		domain = physical_space + 150

		#Endcaps
		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)

		#Water
		sphere3 = Sphere(Pnt(0,0,0),physical_space)
		#PML
		sphere4 = Sphere(Pnt(0,0,0),domain).bc('outer')

		#Delete if it works
		cyl1 = Cylinder(Pnt(0,-2*cyl_length,0),Pnt(0,2*cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-2*physical_space,0),Pnt(0,2*physical_space,0),radius+100)
		cyl3 = Cylinder(Pnt(0,-2*domain,0),Pnt(0,2*domain,0),radius+200).bc('outer')

		box2 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box3 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		middle = cyl1*plane1*plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')
		water = (sphere3 - AuNP).mat('water')
		pmldom = (sphere4 - sphere3).mat('pml').maxh(80)

		geo.Add(AuNP)
		geo.Add(water)
		geo.Add(pmldom)

	else:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		particle = radius + mt_length
		physical_space = length + mt_length + 100
		domain = physical_space + 100

		#Nanoparticle
		cyl1 = Cylinder(Pnt(0,-1,0),Pnt(0,1,0),radius)
		#Microtubules
		cyl2 = Cylinder(Pnt(0,-1,0),Pnt(0,1,0),particle)

		#Nanoparticle endcaps
		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)
		
		#Microtubule endcaps
		sphere3 = Sphere(Pnt(0,-cyl_length,0),particle)
		sphere4 = Sphere(Pnt(0,cyl_length,0),particle)

		#Water 
		sphere5 = Sphere(Pnt(0,0,0),physical_space)
		sphere6 = Sphere(Pnt(0,0,0),domain)

		#Endcap cut-offs
		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		middle = cyl1*plane1*plane2 
		AuNP = (middle + sphere1 + sphere2).mat('gold')

		mt_middle = ((cyl2-AuNP)*plane1*plane2).mat('mt_mid')
		mt_endcap1 = (sphere3- plane1 - sphere1) 
		mt_endcap2 = (sphere4 - plane2 - sphere2)
		mt_endcaps = (mt_endcap1 + mt_endcap2).mat('mt_end')
		
		total = cyl2*plane1*plane2 + sphere3 + sphere4

		water = (sphere5 - total).mat('water')
		pmldom = (sphere6 - sphere5).mat('pml')

		geo.Add(AuNP)
		geo.Add(mt_endcaps)
		geo.Add(mt_middle)
		geo.Add(water)
		geo.Add(pmldom)
		
	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	return mesh