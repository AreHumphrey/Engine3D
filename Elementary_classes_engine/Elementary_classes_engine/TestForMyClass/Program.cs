using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {

           

            Map map = new Map();
            Camera camera = new Camera(new Point(0, 0, 0), new Vector(0, 0, 1), new Vector(0, 0, 1), new Vector(0, 0, 1), 0, 0);
            Point initialPt = new Point(0, 0, 0); Vector dir1 = new Vector(1, 0, 0); Vector dir2 = new Vector(0, 1, 0); Vector dir3 = new Vector(0, 0, 1);
            VectorSpace vectorSpace = new VectorSpace(initialPt, dir1, dir2, dir3);

            
            ParametersSphere sphereParams1 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams2 = new ParametersSphere(1, 1, 1, 7);
            ParametersSphere sphereParams3 = new ParametersSphere(1, 1, 1, 7);
            Sphere sphere1 = new Sphere(new Point(-1, 0, 27), sphereParams1);
            Sphere sphere2 = new Sphere(new Point(6, 0, 40), sphereParams2);
            Sphere sphere3 = new Sphere(new Point(-5, 0, 30), sphereParams3);
            //map.Append(sphere3);
            map.Append(sphere2);
            map.Append(sphere1);
            int numRays = Console.WindowWidth; 
            camera.sendRays(map);
            
            Consoles consoleCanvas = new Consoles(map, camera, vectorSpace);
            consoleCanvas.Draw();
            double fov = 60;
            double drawDistance = 10;
            Spectator spectator = new Spectator(camera, drawDistance);
            Angle a = new Angle(0, 45, 45, 45, 45, 180);
            camera = spectator.Rotate(a);

            Consoles consoleCanvas1 = new Consoles(map, camera, vectorSpace);
            
            consoleCanvas.Draw();
            Console.ReadKey(); 
        }
    }
}

