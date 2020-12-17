from ngsolve import Mesh, pml
from netgen.csg import *

from netgen.geom2d import unit_square

from ngsolve import Draw


#Rename nanoparticle variables if more materials are implemented 

def Nanosphere(particle,mt_length, physical_space, domain, mt_model="radial"):

	# physical_space = particle + mt_length + 150
	# domain = physical_space + 200

	geo = CSGeometry()
	
	if mt_length == 0:

		sphere1 = Sphere(Pnt(0,0,0),particle)
		sphere2 = Sphere(Pnt(0,0,0),physical_space)
		sphere3 = Sphere(Pnt(0,0,0),domain).bc('outer')

		AuNP = sphere1.mat('gold')
		water = (sphere2 - sphere1).mat('water')
		pmldom = (sphere3 - sphere2).mat('pml')

		geo.Add(AuNP)
		geo.Add(water)
		geo.Add(pmldom)

	else:

		if mt_model == "radial":

			sphere1 = Sphere(Pnt(0,0,0),particle)
			sphere2 = Sphere(Pnt(0,0,0),particle+mt_length)
			sphere3 = Sphere(Pnt(0,0,0),physical_space)
			sphere4 = Sphere(Pnt(0,0,0),domain).bc('outer')

			AuNP = sphere1.mat('gold')
			mt = (sphere2 - sphere1).mat('mt_sphere')
			water = (sphere3 - sphere2).mat('water')
			pmldom = (sphere4 - sphere3).mat('pml')

			geo.Add(AuNP)
			geo.Add(mt)
			geo.Add(water)
			geo.Add(pmldom)
		
		else: 

			sphere1 = Sphere(Pnt(0,0,0),particle)
			sphere2 = Sphere(Pnt(0,0,0),physical_space)
			sphere3 = Sphere(Pnt(0,0,0),domain).bc('outer')

			cyl1 = Cylinder(Pnt(0,0,0),Pnt(0,0,1),12.5)
			cyl2 = Cylinder(Pnt(0,0,0),Pnt(0,1,0),12.5)
			cyl3 = Cylinder(Pnt(0,0,0),Pnt(1,0,0),12.5)

			plane1 = Plane(Pnt(0,0,mt_length+particle),Vec(0,0,1))
			plane2 = Plane(Pnt(0,0,-(mt_length+particle)),Vec(0,0,-1))
			plane3 = Plane(Pnt(0,mt_length+particle,0),Vec(0,1,0))
			plane4 = Plane(Pnt(0,-(mt_length+particle),0),Vec(0,-1,0))
			plane5 = Plane(Pnt(mt_length+particle,0,0),Vec(1,0,0))
			plane6 = Plane(Pnt(-(mt_length+particle),0,0),Vec(-1,0,0))

			mt1 = ((plane2*cyl1*plane1) - sphere1).mat("mt_cyl")
			mt2 = ((plane3*cyl2*plane4) - sphere1).mat("mt_cyl")
			mt3 = ((plane5*cyl3*plane6) - sphere1).mat("mt_cyl")

			AuNP = sphere1.mat("gold")
			water = ((((sphere2 - mt1) - mt2) - mt3) - AuNP).mat("water")
			pmldom = (sphere3 - sphere2).mat("pml")

			geo.Add(AuNP)
			geo.Add(mt1)
			geo.Add(mt2)
			geo.Add(mt3)
			geo.Add(water)
			geo.Add(pmldom)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)

	p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
	mesh.SetPML(p,'pml')

	return mesh, particle, particle


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
		pmldom = (cyl3*plane5*plane6 - cyl2*plane3*plane4).mat('pml')

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
		pmldom = ((cyl4-cyl3)*plane5*plane6 - (water + total_body)).mat('pml')

		geo.Add(AuNP)
		geo.Add(mt_endcaps)
		geo.Add(mt_middle)
		geo.Add(water)
		geo.Add(pmldom)
		
	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)

	p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
	mesh.SetPML(p,'pml')

	return mesh, radius, length