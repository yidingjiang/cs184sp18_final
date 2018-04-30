import bpy

LOC = bpy.data.objects["face74.002"].location.copy()
ROT = bpy.data.objects["face74.002"].rotation_euler.copy()
SCALE = bpy.data.objects["face74.002"].scale.copy()
MAT = bpy.data.objects["face74.002"].data.materials[0].copy()

#delete default obj
bpy.ops.object.select_all(action="DESELECT")
bpy.data.objects["face74.002"].select = True
bpy.ops.object.delete()


face_num = 10
for face_num in range(0, 330):
	#import obj
	face_str = 'face'+str(face_num)
	file_loc =  "/home/utkarsh/proj/pbf/cpp/cs184sp18_final/mitsuba/input/" + face_str + ".obj"
	imported_object = bpy.ops.import_scene.obj(filepath=file_loc)
	ob = bpy.context.selected_objects[0]

	#Set obj
	ob.location = LOC
	ob.rotation_euler = ROT
	ob.scale = SCALE

	material_properties = ob.data.materials

	if material_properties:
		material_properties[0] = MAT
	else:
		material_properties.append(MAT)

	#Render scene
	bpy.context.scene.render.filepath = "/home/utkarsh/proj/pbf/cpp/cs184sp18_final/mitsuba/input/img" + str(face_num) + ".png"
	bpy.ops.render.render(write_still=True)

	#Delete obj
	bpy.ops.object.select_all(action="DESELECT")
	bpy.data.objects[face_str].select = True
	bpy.ops.object.delete()


