using System.Diagnostics.Tracing;
using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {
            Map map = new Map();
            Camera camera = new Camera(new Point(0, 0, 0), new Vector(0, 0, 1), new Vector(0, 0, 1), new Vector(0, 0, 1), 80, 100);
            Point initialPt = new Point(0, 0, 0); Vector dir1 = new Vector(1, 0, 0); Vector dir2 = new Vector(0, 1, 0); Vector dir3 = new Vector(0, 0, 1);
            VectorSpace vectorSpace = new VectorSpace(initialPt, dir1, dir2, dir3);
            Event eventSystem = new Event();
            Trigger trigger = new Trigger(eventSystem);

            eventSystem.add("т1");

            eventSystem.handle("т1", (args) => {
                string arg1 = (string)args[0];
                int arg2 = (int)args[1];
                Console.WriteLine($"Событие: {arg1} {arg2}");
            });

            ParametersSphere sphereParams1 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams2 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams3 = new ParametersSphere(1, 1, 1, 7);
            Sphere sphere1 = new Sphere(new Point(-1, 0, 27), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(15, 0, 57), sphereParams2);
            Sphere sphere3 = new Sphere(new Point(-5, 0, 30), sphereParams3);
            //map.Append(sphere3);
            map.Append(sphere2);
            map.Append(sphere1);
            int numRays = Console.WindowWidth; 
            camera.sendRays(map);
            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);
            consoleCanvas.Draw();
            double drawDistance = 10;
            Player pl = new Player(camera, drawDistance, 1, map);
            camera = pl.MoveForward(7);
            consoleCanvas.Draw();
            camera = pl.MoveBack(7);
            camera = pl.Rotate(-3, 0, 0);
            consoleCanvas.Draw();
            
            camera = pl.Rotate(7, 0, 0);
            //camera = pl.MoveForward(20);
            consoleCanvas.Draw();
            trigger.trigger("т1", "т1", 666);
            Console.ReadKey();
            
        }
    }
}

