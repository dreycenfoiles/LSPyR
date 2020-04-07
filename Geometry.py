from netgen.csg import *
from ngsolve import Mesh
from Materials import *

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

		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)

		cyl1 = Cylinder(Pnt(0,-cyl_length,0),Pnt(0,cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-physical_space,0),Pnt(0,physical_space,0),radius+100)
		cyl3 = Cylinder(Pnt(0,-domain,0),Pnt(0,domain,0),radius+200).bc('outer')

		box1 = OrthoBrick(Pnt(-radius,-cyl_length,-radius),Pnt(radius,cyl_length,radius)) 
		box2 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box3 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		middle = (cyl1*box1).mat('gold')
		endcap1 = (sphere1-plane1).mat('gold')
		endcap2 = (sphere2-plane2).mat('gold')

		pmldom = (cyl3*box3).mat('pml').maxh(80)
		hole1 = (pmldom - cyl2*box2)
		water = (pmldom*hole1 - middle - endcap1 - endcap2).mat('water')

		geo.Add(pmldom)
		geo.Add(water)
		geo.Add(middle)
		geo.Add(endcap1)
		geo.Add(endcap2)

	else:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		particle = radius + mt_length
		physical_space = length + mt_length + 100
		domain = physical_space + 100

		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)
		sphere3 = Sphere(Pnt(0,-cyl_length,0),particle)
		sphere4 = Sphere(Pnt(0,cyl_length,0),particle)

		cyl1 = Cylinder(Pnt(0,-cyl_length,0),Pnt(0,cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-cyl_length,0),Pnt(0,cyl_length,0),particle)
		cyl3 = Cylinder(Pnt(0,-physical_space,0),Pnt(0,physical_space,0),particle+100)
		cyl4 = Cylinder(Pnt(0,-domain,0),Pnt(0,domain,0),particle+200).bc('outer')

		box1 = OrthoBrick(Pnt(-radius,-cyl_length,-radius),Pnt(radius,cyl_length,radius)) 
		box2 = OrthoBrick(Pnt(-particle,-cyl_length,-particle),Pnt(particle,cyl_length,particle)) 
		box3 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box4 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))
		
		middle_mt = ((cyl2-cyl1)*box2).mat('mt_mid')
		endcap1_mt = (sphere3-sphere1-plane1).mat('mt_end')
		endcap2_mt = (sphere4-sphere2-plane2).mat('mt_end')

		pmldom = (cyl4*box4).mat('pml').maxh(80)
		hole1 = (pmldom - cyl3*box3)
		water = (pmldom*hole1 - middle_mt - endcap1_mt - endcap2_mt).mat('water')

		middle_rod = (cyl1*plane1*plane2).mat('gold')
		endcap1 = (sphere1-plane1).mat('gold')
		endcap2 = (sphere2-plane2).mat('gold')

		geo.Add(pmldom)
		geo.Add(water)
		geo.Add(middle_mt)
		geo.Add(endcap1_mt)
		geo.Add(endcap2_mt)
		geo.Add(middle_rod)
		geo.Add(endcap1)
		geo.Add(endcap2)


	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	return mesh