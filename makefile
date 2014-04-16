windowsmake: main.o 
	g++ -c -o main.o main.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o controller.o controller/controller.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o bitmap.o bitmap/bitmap.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o graphic.o graphic/graphic.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o object.o object/object.cpp -I "C:\MinGW\freeglut\include"
	g++ -o main.exe main.o controller.o bitmap.o graphic.o object.o -L "C:\MinGW\freeglut\lib" -l freeglut -lopengl32 -lglu32

linuxmake: main.o 
	g++ -c -o main.o main.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o controller.o controller/controller.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o bitmap.o bitmap/bitmap.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o graphic.o graphic/graphic.cpp -I "C:\MinGW\freeglut\include"
	g++ -c -o object.o object/object.cpp -I "C:\MinGW\freeglut\include"
	g++ -o main.exe main.o controller.o bitmap.o graphic.o object.o -L "C:\MinGW\freeglut\lib" -l freeglut -lopengl32 -lglu32

clean:
	rm -f *.o main.exe