using System;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Threading;


namespace Elementary_classes_engine
{
    public class Program
    {
        public class Point
        {
            public double x { get; private set; }
            public double y { get; private set; }
            public double z { get; private set; }

            public Point(double x, double y, double z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }
            public Point Subtract(Vector other) {
return new Point(
this.x - other.x,
this.y - other.y,
this.z - other.z
);
}
            //Перегрузка операторов ( + - * / ) + расстояние
            public static double distance(Point a, Point b)
            {
                return Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
            }
            public static Point operator +(Point a, Point b)
            {
                return new Point(a.x + b.x, a.y + b.y, a.z + b.z);
            }
            public static Point operator -(Point a, Point b)
            {
                return new Point(a.x - b.x, a.y - b.y, a.z - b.z);
            }
            public static Point operator *(Point a, Point b)
            {
                return new Point(a.x * b.x, a.y * b.y, a.z * b.z);
            }
            public static Point operator /(Point a, Point b)
            {
                return new Point(a.x / b.x, a.y / b.y, a.z / b.z);
            }
           

            public static Point operator *(Point a, double element)
            {
                return new Point(a.x * element, a.y * element, a.z * element);
            }
            public static Point operator *(double element, Point b)
            {
                return new Point(b.x * element, b.y * element, b.z * element);
            }
        }
        public class Vector
        {
            public double x = 0, y = 0, z = 0;
            public Point pt;

            public Vector(Point pt)
            {
                this.x = pt.x;
                this.y = pt.y;
                this.z = pt.z;

            }
            public Vector(double x, double y, double z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }

            public static double LenVector(Vector v1)
            {
                return Math.Sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
            }


            //Перегрузка операторов ( + - * / ) + расстояние
            public static Vector operator +(Vector a, Vector b)
            {
                return new Vector(new Point(a.x + b.x, a.y + b.y, a.z + b.z));
            }
            public static Vector operator -(Vector a, Vector b)
            {
                return new Vector(new Point(a.x - b.x, a.y - b.y, a.z - b.z));
            }
            public static Vector operator *(Vector a, Vector b)
            {
                return new Vector(new Point(a.x * b.x, a.y * b.y, a.z * b.z));
            }
            public static Vector operator ^(Vector a, Vector b)
            {
                return new Vector(new Point(a.y * b.z - a.z * b.z, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x));
            }
            public static Vector operator /(Vector a, Vector b)
            {
                return new Vector(new Point(a.x / b.x, a.y / b.y, a.z / b.z));
            }
            public static Vector operator *(double element, Vector b)
            {
                return new Vector(new Point(element * b.x, element * b.y, element * b.z));
            }
            public static Vector operator *(Vector a, double element)
            {
                return new Vector(new Point(element * a.x, element * a.y, element * a.z));
            }
            public static Vector operator /(Vector a, double element)
            {
                return new Vector(new Point(a.x / element, a.y / element, a.z / element));
            }


        }
        public class VectorSpace
        {
            public Point Initpt;
            public Vector dir1;
            public Vector dir2;
            public Vector dir3;

            public VectorSpace(Point InitPt, Vector dir1, Vector dir2, Vector dir3)
            {
                this.Initpt = InitPt;
                this.dir1 = VectorBasis(dir1);
                this.dir2 = VectorBasis(dir2);
                this.dir3 = VectorBasis(dir3);
            }

            public static Vector VectorBasis(Vector vect)
            {
                double vectorsLen = Vector.LenVector(vect);
                return vect / vectorsLen;
            }
        }
        public class Camera
        {
            public double width = 70;
            public double height = 50;
            public double fov;
            public double vFov; //Вертикальный угол отрисовки
            public Point lookAt;
            public Vector lookDir;
            public Point position;
            public Vector rotation;
            public double drawDistance;
            public Camera(Point position, Vector rotation, Point lookAt, Vector lookDir, double fov, double drawDistance)
            {
                this.rotation = rotation;
                this.lookAt = lookAt;
                this.lookDir = lookDir;
                this.fov = fov;
                this.position = position;
               
                this.drawDistance = drawDistance;

            }

            public double DistanceTo(Point point)
            {
                double dx = point.x - position.x;
                double dy = point.y - position.y;
                double dz = point.z - position.z;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }

            public List<Object> nearestObject(Object[] objects)
            {
                List<Object> nearestObjects = new List<Object>();
                double closestDistance = double.PositiveInfinity;
                foreach (Object obj in objects)
                {
                    Ray ray = new Ray(position, lookDir);
                    Point intersection = obj.Intersect(ray);
                    if (intersection != null)
                    {
                        double distance = DistanceTo(intersection);

                        if (distance < closestDistance && distance < drawDistance)
                        {
                            nearestObjects.Add(obj);
                            closestDistance = distance;
                        }
                    }
                }

                return nearestObjects;
            }

            public Vector AngleToDirectionVector(Vector baseDirection, double angle)
            {
                double rad = angle * (Math.PI / 180), cos = Math.Cos(rad), sin = Math.Sin(rad);

                Vector rot = new Vector(baseDirection.x * cos - baseDirection.z * sin, baseDirection.y, baseDirection.x * sin + baseDirection.z * cos);
                return rot;
            }

           
            
            public List<Ray> sendRays(Map map)
            {
                List<Ray> rays = new List<Ray>();

                double rayAngle = fov / width; 
                double vRayAngle = vFov / height; 

                for (int i = 0; i < width; i++)
                {
                    Vector dir = AngleToDirectionVector(lookDir, (i - width / 2) * rayAngle);
                    Ray ray = new Ray(position, dir);
                    rays.Add(ray); 

                    for (int j = 0; j < height; j++)
                    {
                        Vector vDir = AngleToDirectionVector(dir, (j - height / 2) * vRayAngle);
                        Ray vRay = new Ray(position, vDir);
                        rays.Add(vRay); 
                    }
                }

                foreach (Ray ray in rays)
                {
                    Point intersection = ray.Intersect(map);
                    if (intersection != null)
                    {
                        
                    }
                }

                return rays;
            }


        }
       
        public abstract class Object
        {
            public Point position;
            public Vector Rotation;

            public abstract bool Contains(Point Pt);

            public  virtual Point Intersect(Ray ray)
            {
                return null;
            }

            public virtual Point NearestPoint(Point[] points)
            {
                Point nearestPoint = null; 
                double closestDistance = double.PositiveInfinity; 

                foreach (Point point in points)
                {
                    double distance = DistanceTo(point);
                    if (distance < closestDistance)
                    {
                        nearestPoint = point;
                        closestDistance = distance;
                    }
                }

                return nearestPoint;
            }

            public double DistanceTo(Point point)
            {
                double dx = point.x - position.x;
                double dy = point.y - position.y;
                double dz = point.z - position.z;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }



        }

        public class Plane : Object // Плоскость наследуется от объекта 
        {
            public Parameters parameters;

            public Plane(Point position, Vector rotation)
            {
                this.position = position;
                this.Rotation = rotation;
                parameters = new Parameters();
            }

            public override bool Contains(Point point)
            {
                double result = parameters.aCoefficient * point.x + parameters.bCoefficient * point.y + parameters.cCoefficient * point.z;

                return result == 0;
            }

            public override Point Intersect(Ray ray)
            {
                double denominator = parameters.aCoefficient * ray.direction.x + parameters.bCoefficient * ray.direction.y + parameters.cCoefficient * ray.direction.z;
                if (denominator == 0)
                {
                    if (Contains(ray.position))
                    {
                        return ray.position;
                    }
                    else
                    {
                        return null;
                    }
                }
                else
                {
                    double numerator = parameters.aCoefficient * (ray.position.x - position.x) + parameters.bCoefficient * (ray.position.y - position.y) + parameters.cCoefficient * (ray.position.z - position.z);

                    double distance = numerator / denominator;

                    return new Point(ray.position.x + ray.direction.x * distance, ray.position.y + ray.direction.y * distance, ray.position.z + ray.direction.z * distance);
                }
            }
        }
        public class Parameters
        {
            public double aCoefficient;
            public double bCoefficient;
            public double cCoefficient;

            public void Rotate(double angle) 
            {
                double cosA = Math.Cos(angle);
                double sinA = Math.Sin(angle);

                double newACoefficient = aCoefficient * cosA - bCoefficient * sinA;
                double newBCoefficient = bCoefficient * cosA + aCoefficient * sinA;

                aCoefficient = newACoefficient;
                bCoefficient = newBCoefficient;
            }

            public void Move(double x, double y, double z) // Преобразовать коэффициенты уравнения при помощи перемещения
            {
                aCoefficient += x;
                bCoefficient += y;
                cCoefficient += z;
            }

            public void Scale(double xFactor, double yFactor, double zFactor)  // Преобразовать коэффициенты уравнения при помощи масштабирования
            {
                aCoefficient *= xFactor;
                bCoefficient *= yFactor;
                cCoefficient *= zFactor;
            }
        }
        public class ParametersPlane : Parameters
        {
            public double a;
            public double b;
            public double c;
            public double d;   
            public ParametersPlane()
            {
                a = 0;
                b = 0;
                c = 1;
                d = 0;
            }
            
        }
        public class BoundedPlane : Plane //  Главная плоскость
        {
            private double dx;
            private double dy;
            private double dz;

            public BoundedPlane(Point position, Vector rotation, double dx, double dy, double dz) : base(position, rotation)
            {
                this.dx = dx;
                this.dy = dy;
                this.dz = dz;

            }
            public void InBoundaries(Point Pt)  // Проверка координат точки на соответствие границам плоскости
            {
           
            }

        }


        public class Sphere : Object
        {
            private ParametersSphere parameters;


            public Sphere(Point Center, ParametersSphere parameters)
            {
                position = Center;

                this.parameters = parameters;
            }

            public override bool Contains(Point point)
            {
                return Math.Pow(point.x - position.x, 2) + Math.Pow(point.y - position.y, 2) + Math.Pow(point.z - position.z, 2) == Math.Pow(parameters.radius, 2);
            }

            public override Point Intersect(Ray ray)
            {
                double a = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y + ray.direction.z * ray.direction.z,
                b = 2 * (ray.direction.x * (ray.position.x - position.x) + ray.direction.y * (ray.position.y - position.y) +  ray.direction.z * (ray.position.z - position.z)),
                c = (ray.position.x - position.x) * (ray.position.x - position.x) + (ray.position.y - position.y) * (ray.position.y - position.y) + (ray.position.z - position.z) * (ray.position.z - position.z) - parameters.radius * parameters.radius;

                double discriminant = b * b - 4 * a * c;

                if (discriminant < 0)
                {
                    return null;
                }

                double t1 = (-b + Math.Sqrt(discriminant)) / (2 * a), t2 = (-b - Math.Sqrt(discriminant)) / (2 * a);

                if (Math.Abs(t1) < double.Epsilon && Math.Abs(t2) < double.Epsilon)
                {
                    return null;
                }

                double t = Math.Abs(t1) > Math.Abs(t2) ? t2 : t1;

                return new Point(ray.position.x + t * ray.direction.x, ray.position.y + t * ray.direction.y, ray.position.z + t * ray.direction.z);

            }



        }
        public class ParametersSphere : Parameters
        {
            public double radius;

            public ParametersSphere(double x, double y, double z, double radius) 
            {
                this.radius = radius;
                this.aCoefficient = x;
                this.bCoefficient = y;
                this.cCoefficient = z;

            }
        }
        public class Map
        {
            public Object[] arrObj;
            public void Append(Object obj)
            {
                if (arrObj == null)
                {
                    arrObj = new Object[] { obj };
                }
                else
                {
                    Object[] newArrObjects = new Object[arrObj.Length + 1];
                    for (int i = 0; i < arrObj.Length; i++)
                    {
                        newArrObjects[i] = arrObj[i];
                    }
                    newArrObjects[newArrObjects.Length - 1] = obj;
                    arrObj = newArrObjects;
                }
            }

        }
        public class Ray
        {
            public Point position;
            public Vector direction;

            public Ray(Point position, Vector direction)
            {
                this.position = position;
                this.direction = direction;
            }

            public Point Intersect(Map map)
            {
                Object nearestObject = null;
                double nearestDistance = double.PositiveInfinity;
                Point nearestPoint = null;

                foreach (Object obj in map.arrObj)
                {
                    Point point = obj.Intersect(this); // передаем объекту текущий луч
                    if (point != null)
                    {
                        double distance = obj.DistanceTo(position);
                        if (distance < nearestDistance)
                        {
                            nearestObject = obj;
                            nearestDistance = distance;
                            nearestPoint = point;
                        }
                    }
                }

                return nearestPoint;
            }


        }
        
        public class Canvas
        {
            public Map map;
            public Camera camera;
            public VectorSpace vSpace;
            public Canvas(Map map, Camera camera, VectorSpace vectorSpace)
            {
                this.map = map;
                this.camera = camera;
                this.vSpace = vectorSpace;
            }
            public void draw()
            {

            }  
        }

        public class Consoles : Canvas
        {
            private ConsoleColor[] gradient = new ConsoleColor[]
            {
                ConsoleColor.Red, ConsoleColor.Green, ConsoleColor.Blue, ConsoleColor.Yellow, ConsoleColor.DarkYellow, ConsoleColor.DarkCyan

            };
            public Consoles(Map map, Camera camera, VectorSpace vSpace) : base(map, camera, vSpace) { }

            private string dropline = "@#*+=-^:. ";
  
            public void Draw()
            {
                int screenWidth = 100;
                int screenHeight = 50;

                for (int y = 0; y < screenHeight; y++)
                {
                    for (int x = 0; x < screenWidth ; x++)
                    {

                        Vector dir = vSpace.dir1 * ((double)x / screenWidth - 0.5) + vSpace.dir2 * ((double)y / screenHeight - 0.5) + vSpace.dir3;

                        List<Object> nearestObjects = camera.nearestObject(map.arrObj);
                        foreach (var obj in nearestObjects)
                        {
                            if (obj != null)
                            {
                                Ray ray = new Ray(camera.position, dir);
                                Point intersection = obj.Intersect(ray);

                                if (intersection != null)
                                {
                                    double distance = camera.DistanceTo(intersection);
                                    double gradientIndex = (distance / 100.0) * dropline.Length;
                                    for (int i = 1; i < obj.DistanceTo(camera.position); i++)
                                    {
                                        gradientIndex += 0.05; // изменение на 0.5% за каждый пиксель
                                    }
                                    gradientIndex = Math.Min(gradientIndex, dropline.Length - 1);
                                    Console.Write(dropline[(int)Math.Round(gradientIndex)]);
                                }
                                else
                                {
                                    Console.Write(" ");
                                }
                            }
                            else
                            {
                                Console.Write(" ");
                            }
                        }


                    }

                    Console.WriteLine();
                }



            }

        }
        public class Angle : Point //Угол поворота 
        {
            public double yz { get; set; } 
            public double xz { get; set; }
            public double xy { get; set; }
            public Angle(double x, double y, double z, double yz, double xz, double xy) : base(x, y, z)
            {
                this.yz = ToRange(yz);
                this.xz = ToRange(xz);
                this.xy = ToRange(xy);
            }
            public double ToRange(double angle)
            {
                angle %= 360;
                if (angle < 0)
                {
                    angle += 360;
                }
                return angle;
            }
            public static Angle operator +(Angle a, Angle b)
            {
                double newYZ = a.yz + b.yz;
                double newXZ = a.xz + b.xz;
                double newXY = a.xy + b.xy;
                return new Angle(0, 0, 0, newYZ, newXZ, newXY);
            }

            public static Angle operator -(Angle a)
            {
                double newYz = -a.yz;
                double newXz = -a.xz;
                double newXy = -a.xy;
                return new Angle(0, 0, 0, newYz, newXz, newXy);
            }

            public static Vector operator *(Point a, Angle angle)
            {
                double yzRad = angle.yz * Math.PI / 180;
                double xzRad = angle.xz * Math.PI / 180;
                double xyRad = angle.xy * Math.PI / 180;

                double[,] matrix = new double[3, 3];
                matrix[0, 0] = Math.Cos(xzRad) * Math.Cos(yzRad);
                matrix[0, 1] = Math.Cos(xzRad) * Math.Sin(yzRad);
                matrix[0, 2] = -Math.Sin(xzRad);

                matrix[1, 0] = -Math.Cos(xyRad) * Math.Sin(yzRad) + Math.Sin(xyRad) * Math.Sin(xzRad) * Math.Cos(yzRad);
                matrix[1, 1] = Math.Cos(xyRad) * Math.Cos(yzRad) + Math.Sin(xyRad) * Math.Sin(xzRad) * Math.Sin(yzRad);
                matrix[1, 2] = Math.Sin(xyRad) * Math.Cos(xzRad);

                matrix[2, 0] = Math.Sin(xyRad) * Math.Sin(yzRad) + Math.Cos(xyRad) * Math.Sin(xzRad) * Math.Cos(yzRad);
                matrix[2, 1] = -Math.Sin(xyRad) * Math.Cos(yzRad) + Math.Cos(xyRad) * Math.Sin(xzRad) * Math.Sin(yzRad);
                matrix[2, 2] = Math.Cos(xyRad) * Math.Cos(xzRad);

                double newX = matrix[0, 0] * a.x + matrix[0, 1] * a.y + matrix[0, 2] * a.z;
                double newY = matrix[1, 0] * a.x + matrix[1, 1] * a.y + matrix[1, 2] * a.z;
                double newZ = matrix[2, 0] * a.x + matrix[2, 1] * a.y + matrix[2, 2] * a.z;

                Point t = new Point(newX, newY, newZ);
                return new Vector(t);

            }


        }

        

        public class Player
        {
            public Camera camera;
            public Point position;
            public Vector lookAt;
            public double collisionRadius;
            public Map map;

            public Player(Point position, double fov, double drawDistance, double collisionRadius, Map map)
            {
                this.camera = new Camera(position, new Vector(0, 0, 0), null, null, fov, drawDistance);
                this.position = position;
                this.lookAt = new Vector(0, 0, 1); // default look direction is towards +Z axis
                this.collisionRadius = collisionRadius;
                this.map = map;
            }

            public void MoveForward(double distance)
            {
                double dx = distance * lookAt.x;
                double dz = distance * lookAt.z;
                Point newPosition = new Point(position.x + dx, position.y, position.z + dz);
                if (!CheckCollision(newPosition))
                {
                    position = newPosition;
                    camera.position = newPosition;
                }
            }

            public void MoveBack(double distance)
            {
                double dx = distance * lookAt.x;
                double dz = distance * lookAt.z;
                Point newPosition = new Point(position.x - dx, position.y, position.z - dz);
                if (!CheckCollision(newPosition))
                {
                    position = newPosition;
                    camera.position = newPosition;
                }
            }

            public void StrafeLeft(double distance)
            {
                double dx = distance * lookAt.z;
                double dz = -distance * lookAt.x;
                Point newPosition = new Point(position.x + dx, position.y, position.z + dz);
                if (!CheckCollision(newPosition))
                {
                    position = newPosition;
                    camera.position = newPosition;
                }
            }

            public void StrafeRight(double distance)
            {
                StrafeLeft(-distance);
            }

            public void Rotate(Angle angle)
            {
                Vector newRotation = new Vector(camera.rotation.x + angle.xz, camera.rotation.y + angle.yz, camera.rotation.z + angle.xy);
                camera.rotation = newRotation;
                lookAt = new Vector(Math.Sin(newRotation.x * Math.PI / 180), Math.Tan(newRotation.y * Math.PI / 180), Math.Cos(newRotation.x * Math.PI / 180) * Math.Cos(newRotation.y * Math.PI / 180));
            }

            private bool CheckCollision(Point newPosition)
            { 
                foreach (Object obj in map.arrObj) 
                {
                    if (obj.Contains(newPosition))
                    { 
                        return true; 
                    }
                }
                return false;
            }
        }
        public class FreeCamera : Spectator 
        {
            public FreeCamera(Point position, Vector rotation, double fov, double drawDistance) : base(position, rotation, fov, drawDistance) { }

        }
        public class Spectator
        {
            public Camera camera;
            public Point position;
            public Vector lookAt;
            //public Spectator( Camera camera, Point position, Vector lookAt)
            //{
            //    this.camera = camera;
            //    this.position = position;
            //    this.lookAt = lookAt;
            //}
            public Spectator(Point position, Vector rotation, double fov, double drawDistance)
            {
                this.camera = new Camera(position, rotation, null, null, fov, drawDistance);
            }

            public Spectator(Point position, double fov, double drawDistance)
                : this(position, new Vector(0, 0, 0), fov, drawDistance) { }

            public void MoveBack(double distance) 
            {
                double dx = distance * camera.lookDir.x;
                double dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x - dx, camera.position.y, camera.position.z - dz);
            }
            public void MaveForward(double distance)
            {
                //position = position + camera.position * distance;
                //camera.position = position;
                double dx = distance * camera.lookDir.x; 
                double dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);
            }
            public void StrafeLeft(double distance)
            {
                double dx = distance * camera.lookDir.z;
                double dz = -distance * camera.lookDir.x;
                camera.position = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);
            }
            public void StrafeRight(double distance)
            {
                StrafeLeft(-distance);
            }
            public void Rotate(Angle angle)
            {
                Vector newRotation = new Vector(camera.rotation.x + angle.xz, camera.rotation.y + angle.yz, camera.rotation.z + angle.xy);
                camera.rotation = newRotation;
                camera.lookDir = new Vector(Math.Sin(newRotation.x * Math.PI / 180), Math.Tan(newRotation.y * Math.PI / 180), Math.Cos(newRotation.x * Math.PI / 180) * Math.Cos(newRotation.y * Math.PI / 180));
            }


        }
       

        





    }


}

