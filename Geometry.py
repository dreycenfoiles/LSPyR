from netgen.csg import *
from ngsolve import Mesh


#Rename nanoparticle variables if more materials are implemented 

def Nanosphere(particle,mt_length):

	physical_space = particle + mt_length + 150
	domain = physical_space + 100

	geo = CSGeometry()
	
	if mt_length == 0:

		sphere1 = Sphere(Pnt(0,0,0),particle)
		sphere2 = Sphere(Pnt(0,0,0),physical_space)
		sphere3 = Sphere(Pnt(0,0,0),domain).bc('outer')

		AuNP = sphere1.mat('gold')
		water = (sphere2 - sphere1).mat('water')
		pml = (sphere3 - sphere2).mat('pml').maxh(80)

		geo.Add(AuNP)
		geo.Add(water)
		geo.Add(pml)

	else:

		sphere1 = Sphere(Pnt(0,0,0),particle)
		sphere2 = Sphere(Pnt(0,0,0),particle+mt_length)
		sphere3 = Sphere(Pnt(0,0,0),physical_space)
		sphere4 = Sphere(Pnt(0,0,0),domain).bc('outer')

		AuNP = sphere1.mat('gold')
		mt = (sphere2 - sphere1).mat('mt_sphere')
		water = (sphere3 - sphere2).mat('water')
		pml = (sphere4 - sphere3).mat('pml').maxh(80)

		geo.Add(AuNP)
		geo.Add(mt)
		geo.Add(water)
		geo.Add(pml)

	ngmesh = geo.GenerateMesh(maxh=40)

	mesh = Mesh(ngmesh)

	return mesh


def Nanorod(aeff,ratio,mt_length):

	geo = CSGeometry()

	if ratio == 1:
		return Nanosphere(aeff,mt_length)

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

	
		cyl1 = Cylinder(Pnt(0,-2*cyl_length,0),Pnt(0,2*cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-2*physical_space,0),Pnt(0,2*physical_space,0),radius+100)
		cyl3 = Cylinder(Pnt(0,-2*domain,0),Pnt(0,2*domain,0),radius+200).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		plane3 = Plane(Pnt(0,-physical_space,0),Vec(0,-1,0))
		plane4 = Plane(Pnt(0,physical_space,0),Vec(0,1,0))

		plane5 = Plane(Pnt(0,-domain,0),Vec(0,-1,0))
		plane6 = Plane(Pnt(0,domain,0),Vec(0,1,0))

		middle = cyl1*plane1*plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')
		water = (cyl2*plane3*plane4 - AuNP).mat('water')
		pmldom = (cyl3*plane5*plane6 - cyl2*plane3*plane4).mat('pml').maxh(80)

		geo.Add(AuNP)
		geo.Add(water)
		geo.Add(pmldom)

	else:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		particle = radius + mt_length
		physical_space = length + mt_length + 100
		domain = physical_space + 150

		

		#Nanoparticle
		cyl1 = Cylinder(Pnt(0,-1,0),Pnt(0,1,0),radius)
		#Microtubules
		cyl2 = Cylinder(Pnt(0,-1,0),Pnt(0,1,0),particle)
		#Water
		cyl3 = Cylinder(Pnt(0,-1,0),Pnt(0,1,0),particle+100)
		#PML
		cyl4 = Cylinder(Pnt(0,-2,0),Pnt(0,2,0),particle+150).bc('outer')

		#Nanoparticle endcaps
		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)
		
		#Microtubule endcaps
		sphere3 = Sphere(Pnt(0,-cyl_length,0),particle)
		sphere4 = Sphere(Pnt(0,cyl_length,0),particle)


		#Endcap cut-offs
		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		plane3 = Plane(Pnt(0,-physical_space,0),Vec(0,-1,0))
		plane4 = Plane(Pnt(0,physical_space,0),Vec(0,1,0))

		plane5 = Plane(Pnt(0,-domain,0),Vec(0,-1,0))
		plane6 = Plane(Pnt(0,domain,0),Vec(0,1,0))



		middle = cyl1*plane1*plane2 
		AuNP = (middle + sphere1 + sphere2).mat('gold')

		mt_middle = ((cyl2-cyl1)*plane1*plane2).mat('mt_mid')
		endcap1_mt = (sphere3-plane1)-sphere1
		endcap2_mt = (sphere4-plane2)-sphere2
		mt_endcaps = (endcap1_mt + endcap2_mt).mat('mt_end')
		total_body = AuNP + endcap1_mt + endcap2_mt + mt_middle

		water = (cyl3*plane3*plane4 - total_body).mat('water')
		pmldom = ((cyl4-cyl3)*plane5*plane6 - (water + total_body)).mat('pml').maxh(60)

		geo.Add(AuNP)
		geo.Add(mt_endcaps)
		geo.Add(mt_middle)
		geo.Add(water)
		geo.Add(pmldom)
		
	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	return mesh