using Microsoft.Win32.SafeHandles;
using System;
using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {
            
            Map map = new Map();
            Camera camera = new Camera(new Point(0, 0, 0), new Vector(0, 0, 1), new Vector(0, 0, 1), new Vector(0, 0, 1), 0, 100);
            Point initialPt = new Point(0, 0, 0); Vector dir1 = new Vector(1, 0, 0); Vector dir2 = new Vector(0, 1, 0); Vector dir3 = new Vector(0, 0, 1);
            VectorSpace vectorSpace = new VectorSpace(initialPt, dir1, dir2, dir3);
            Event eventSystem = new Event();
            Trigger trigger = new Trigger(eventSystem);
            ParametersSphere sphereParams1 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams2 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams3 = new ParametersSphere(1, 1, 1, 7);
            Sphere sphere1 = new Sphere(new Point(-1, -5, 27), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(15, 0, 47), sphereParams2);
            Sphere sphere3 = new Sphere(new Point(-5, 0, 40), sphereParams3);
           
            Angle angle = new Angle(5, 20, 10, 5, 20, 10);

            ParametersParallelepiped parParams = new ParametersParallelepiped(0, 0, 0, 5, 5, 5);
            Parallelepiped parallelepiped = new Parallelepiped(new Point(5, 0, 20), parParams);

            map.Append(parallelepiped);
           
            //map.Append(sphere3);
            //map.Append(sphere1);
            //map.Append(sphere2);
            
            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);
            consoleCanvas.Draw();
            double drawDistance = 10;
            Player pl = new Player(camera, drawDistance, 1, map);
            eventSystem.add("MoveForward");
            eventSystem.add("MoveBack");
            eventSystem.add("StrafeLeft");
            eventSystem.add("StrafeRight");
            eventSystem.add("Rotate");
            eventSystem.add("MoveUp");
            eventSystem.add("MoveDown");
            while (true)
            {
                if (Console.KeyAvailable)
                {
                    ConsoleKeyInfo keyInfo = Console.ReadKey(true);
                    switch (keyInfo.Key)
                    {
                        case ConsoleKey.W:
                            camera = pl.MoveForward(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.S:
                            camera = pl.MoveBack(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.A:
                            camera = pl.StrafeLeft(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.D:
                         
                            camera = pl.StrafeRight(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.Q:
                            camera = pl.Rotate(-2, 0, 0);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.E:
                            camera = pl.Rotate(2, 0, 0);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.Z:
                           camera = pl.MoveUp(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;
                        case ConsoleKey.X:  
                            camera = pl.MoveDown(0.5);
                            Console.Clear();
                            consoleCanvas.Draw();
                            System.Threading.Thread.Sleep(20);
                            break;

                        

                    }
                }
            }
        }
    }
}