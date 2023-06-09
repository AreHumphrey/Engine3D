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
            ParametersSphere sphereParams3 = new ParametersSphere(1, 1, 1, 5);
            ParametersSphere sphereParams4 = new ParametersSphere(1, 1, 1, 7);
            Sphere sphere1 = new Sphere(new Point(-1, -7, 40), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(15, -9, 47), sphereParams2);
            Sphere sphere3 = new Sphere(new Point(15, 5, 25), sphereParams3);
            Sphere sphere4 = new Sphere(new Point(2, 0, 35), sphereParams4);
            Angle angle = new Angle(5, 20, 10, 5, 20, 10);
            ParametersParallelepiped parParams1 = new ParametersParallelepiped(0, 0, 0, 2, 10, 10);
            Parallelepiped parallelepiped1 = new Parallelepiped(new Point(0, -2, 20), parParams1);
            ParametersParallelepiped parParams2 = new ParametersParallelepiped(-5, -2, 0, 8, 8, 8);
            Parallelepiped parallelepiped2 = new Parallelepiped(new Point(4, 4, 40), parParams2);
            ParametersParallelepiped parParams3 = new ParametersParallelepiped(-4, 2, 0, 12, 12, 12);
            Parallelepiped parallelepiped3 = new Parallelepiped(new Point(-6, 7, 28), parParams2);
            map.Append(parallelepiped1);
            //map.Append(parallelepiped2);
            //map.Append(parallelepiped3);
            //map.Append(sphere4);
            //map.Append(sphere3);
            //map.Append(sphere2);
            //map.Append(sphere1);
            
            

            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);
            consoleCanvas.Draw();
            double drawDistance = 10;
            Player pl = new Player(camera, drawDistance, 1, map);
            eventSystem.Append("MoveForward");
            eventSystem.Handle("MoveForward", (args) =>
            {
                camera = pl.MoveForward((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });
            eventSystem.Append("MoveBack");
            eventSystem.Handle("MoveBack", (args) =>
            {
                camera = pl.MoveBack((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });
            eventSystem.Append("StrafeLeft");
            eventSystem.Handle("StrafeLeft", (args) =>
            {
                camera = pl.StrafeLeft((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });
            eventSystem.Append("StrafeRight");
            eventSystem.Handle("StrafeRight", (args) =>
            {
                camera = pl.StrafeRight((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });
            eventSystem.Append("Rotate");
            eventSystem.Append("MoveUp");
            eventSystem.Handle("MoveUp", (args) =>
            {
                camera = pl.MoveUp((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });
            eventSystem.Append("MoveDown");
            eventSystem.Handle("MoveDown", (args) =>
            {
                camera = pl.MoveDown((double)args[0]);
                Console.Clear();
                consoleCanvas.Draw();
                System.Threading.Thread.Sleep(20);
            });

            while (true)
            {
                if (Console.KeyAvailable)
                {
                    ConsoleKeyInfo keyInfo = Console.ReadKey(true);
                    switch (keyInfo.Key)
                    {
                        case ConsoleKey.W:
                            trigger.Triggers("MoveForward", 0.5);
                            break;
                        case ConsoleKey.S:
                            trigger.Triggers("MoveBack", 0.5);
                            break;
                        case ConsoleKey.A:
                            trigger.Triggers("StrafeLeft", 0.5);
                            break;
                        case ConsoleKey.D:
                            trigger.Triggers("StrafeRight", 0.5);
                            break;
                        case ConsoleKey.Z:
                            trigger.Triggers("MoveUp", 0.5);
                            break;
                        case ConsoleKey.X:
                            trigger.Triggers("MoveDown", 0.5);
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
                    }
                }
            }
 
        }
    }
}
