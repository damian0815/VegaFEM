all:
	echo "Compiling all libraries and driver..."
	make -C utilities/interactiveDeformableSimulator 
	
	echo "Compiling volumetric mesh utilities..."
	make -C utilities/volumetricMeshUtilities


clean:
	make -C utilities/interactiveDeformableSimulator clean
	make -C utilities/volumetricMeshUtilities deepclean

