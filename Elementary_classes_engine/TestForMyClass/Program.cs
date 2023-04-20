
using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {

            Console.SetWindowSize(100, 40);

           
            Map map = new Map();
            Camera camera = new Camera(new Point(0, 0, 0), new Vector(0, 0, 1), new Point(0, 0, 1), new Vector(0, 0, 1), 0, 50);
            Point initialPt = new Point(0, 0, 0); Vector dir1 = new Vector(1, 0, 0); Vector dir2 = new Vector(0, 1, 0); Vector dir3 = new Vector(0, 0, 1);
            VectorSpace vectorSpace = new VectorSpace(initialPt, dir1, dir2, dir3);

            
            ParametersSphere sphereParams1 = new ParametersSphere(1, 1, 1, 4);
            ParametersSphere sphereParams2 = new ParametersSphere(1, 1, 1, 2);
            Sphere sphere1 = new Sphere(new Point(0, 0, 12), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(0, 0, 20), sphereParams1);
            map.Append(sphere1);
            map.Append(sphere2);

            int numRays = Console.WindowWidth;
            int rays = Console.WindowHeight;
            camera.SendRays(map, numRays, rays);

            
            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);

            
            consoleCanvas.Draw();

            Console.ReadKey();




        }
    }
}

