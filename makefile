CC = g++
INC_DIR = "C:\MinGW\freeglut\include"
VPATH = object/ bitmap/ graphic/ tm/

OBJPATH = obj/
OBJFLAGS = -c -o


bygg: 
	make -s $(OBJPATH)main.o
	make -s $(OBJPATH)controller.o
	make -s $(OBJPATH)bitmap.o
	make -s $(OBJPATH)graphic.o
	make -s $(OBJPATH)object.o
	make -s $(OBJPATH)object.o
	make -s $(OBJPATH)objectSubdivide.o
	make -s $(OBJPATH)objectTruncate.o
	make -s $(OBJPATH)objectExpand.o
	make -s $(OBJPATH)objSimple.o
	make -s $(OBJPATH)objCube.o
	make -s $(OBJPATH)objDodecahedron.o
	make -s $(OBJPATH)objTetrahedron.o
	make -s $(OBJPATH)objTruncatedIcosahedron.o
	$(CC) -o geo.exe $(OBJPATH)*.o -L "C:\MinGW\freeglut\lib" -l freeglut -lopengl32 -lglu32

$(OBJPATH)main.o: main.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)main.o main.cpp -I $(INC_DIR)

$(OBJPATH)controller.o: controller/controller.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)controller.o controller/controller.cpp -I $(INC_DIR)

$(OBJPATH)bitmap.o: bitmap/bitmap.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)bitmap.o bitmap/bitmap.cpp -I $(INC_DIR)

$(OBJPATH)graphic.o: graphic/graphic.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)graphic.o graphic/graphic.cpp -I $(INC_DIR)

$(OBJPATH)object.o: object/object.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)object.o object/object.cpp -I $(INC_DIR)

$(OBJPATH)objectSubdivide.o: object/objectSubdivide.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objectSubdivide.o object/objectSubdivide.cpp -I $(INC_DIR)

$(OBJPATH)objectTruncate.o: object/objectTruncate.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objectTruncate.o object/objectTruncate.cpp -I $(INC_DIR)

$(OBJPATH)objectExpand.o: object/objectExpand.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objectExpand.o object/objectExpand.cpp -I $(INC_DIR)	

$(OBJPATH)objSimple.o: object/objSimple.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objSimple.o object/objSimple.cpp -I $(INC_DIR)

$(OBJPATH)objCube.o: object/objCube.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objCube.o object/objCube.cpp -I $(INC_DIR)

$(OBJPATH)objDodecahedron.o: object/objDodecahedron.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objDodecahedron.o object/objDodecahedron.cpp -I $(INC_DIR)

$(OBJPATH)objTetrahedron.o: object/objTetrahedron.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objTetrahedron.o object/objTetrahedron.cpp -I $(INC_DIR)

$(OBJPATH)objTruncatedIcosahedron.o: object/objTruncatedIcosahedron.cpp
	$(CC) $(OBJFLAGS) $(OBJPATH)objTruncatedIcosahedron.o object/objTruncatedIcosahedron.cpp -I $(INC_DIR)

object2:
	make $(OBJPATH)*.o

objects: main.cpp

	

clean: 
	rm -f $(OBJPATH)*.o geo.exe