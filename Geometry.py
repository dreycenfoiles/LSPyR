from netgen.csg import *
from ngsolve import *

def Nanosphere():

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

		cyl1 = Cylinder(Pnt(0,-2*cyl_length,0),Pnt(0,2*cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-2*physical_space,0),Pnt(0,2*physical_space,0),radius+100)
		cyl3 = Cylinder(Pnt(0,-2*domain,0),Pnt(0,2*domain,0),radius+200).bc('outer')

		box2 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box3 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		middle = cyl1*plane1*plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')
		water = (cyl2*box2 - AuNP).mat('water')
		pmldom = (cyl3*box3 - cyl2*box2).mat('pml').maxh(80)

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
		
		middle = cyl1*plane1*plane2
		endcap1 = sphere1-plane1
		endcap2 = sphere2-plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')

		middle_mt = ((cyl2-cyl1)*plane1*plane2).mat('mt_mid')
		endcap1_mt = (sphere3-plane1)-sphere1
		endcap2_mt = (sphere4-plane2)-sphere2
		endcaps_mt = (endcap1_mt + endcap2_mt).mat('mt_end')
		total_body = AuNP + middle_mt + endcaps_mt

		water = (cyl3*box3 - total_body).mat('water')
		pmldom = (cyl4*box4 - cyl3*box3).mat('pml').maxh(80)

		geo.Add(AuNP)
		geo.Add(middle_mt)
		geo.Add(endcaps_mt)
		geo.Add(water)
		geo.Add(pmldom)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	return mesh