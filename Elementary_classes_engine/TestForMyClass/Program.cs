﻿using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {

            Console.SetWindowSize(100, 50);

            Map map = new Map();
            Camera camera = new Camera(new Point(0, 0, 0), new Vector(0, 0, 1), new Point(0, 0, 1), new Vector(0, 0, 1), 0, 80);
            Point initialPt = new Point(0, 0, 0); Vector dir1 = new Vector(1, 0, 0); Vector dir2 = new Vector(0, 1, 0); Vector dir3 = new Vector(0, 0, 1);
            VectorSpace vectorSpace = new VectorSpace(initialPt, dir1, dir2, dir3);

            
            ParametersSphere sphereParams1 = new ParametersSphere(1, 1, 1, 5);
            ParametersSphere sphereParams2 = new ParametersSphere(1, 1, 1, 6);
            Sphere sphere1 = new Sphere(new Point(0, 0, 36), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(1, 0, 100), sphereParams2);

            map.Append(sphere2);
            map.Append(sphere1);
            
            
            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);
           
            consoleCanvas.Draw();

            Console.ReadKey(); 
        }
    }
}

