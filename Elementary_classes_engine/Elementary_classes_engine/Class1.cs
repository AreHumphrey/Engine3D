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
           

            public static Point operator *(Point a, double Element)
            {
                return new Point(a.x * Element, a.y * Element, a.z * Element);
            }
            public static Point operator *(double Element, Point b)
            {
                return new Point(b.x * Element, b.y * Element, b.z * Element);
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
            public static Vector operator *(double Element, Vector b)
            {
                return new Vector(new Point(Element * b.x, Element * b.y, Element * b.z));
            }
            public static Vector operator *(Vector a, double Element)
            {
                return new Vector(new Point(Element * a.x, Element * a.y, Element * a.z));
            }
            public static Vector operator /(Vector a, double element)
            {
                return new Vector(new Point(a.x / element, a.y / element, a.z / element));
            }


        }
        public class VectorSpace
        {
            public Point Initpt;
            public Vector Dir1;
            public Vector Dir2;
            public Vector Dir3;

            public VectorSpace(Point InitPt, Vector Dir1, Vector Dir2, Vector Dir3)
            {
                this.Initpt = InitPt;
                this.Dir1 = VectorBasis(Dir1);
                this.Dir2 = VectorBasis(Dir2);
                this.Dir3 = VectorBasis(Dir3);
            }

            public static Vector VectorBasis(Vector vect)
            {
                double vectorsLen = Vector.LenVector(vect);
                return vect / vectorsLen;
            }
        }
        public class Camera
        {
            public double Width = 70;
            public double Height = 50;
            public double Fov;
            public double VFov; //Вертикальный угол отрисовки
            public Point Look_at;
            public Vector Look_dir;
            public Point position;
            public Vector rotation;
            public double DrawDistance;
            public Camera(Point position, Vector rotation, Point Look_at, Vector Look_dir, double Fov, double DrawDistance)
            {
                this.rotation = rotation;
                this.Look_at = Look_at;
                this.Look_dir = Look_dir;
                this.Fov = Fov;
                this.position = position;
               
                this.DrawDistance = DrawDistance;

            }

            public double DistanceTo(Point point)
            {
                double dx = point.x - position.x;
                double dy = point.y - position.y;
                double dz = point.z - position.z;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }

            public Object NearestObject(Object[] objects)
            {
                Object nearestObject = null;
                double closestDistance = double.PositiveInfinity;

                foreach (Object obj in objects)
                {
                    Ray ray = new Ray(position, Look_dir);
                    Point intersection = obj.Intersect(ray);
                    if (intersection != null)
                    {
                        double distance = DistanceTo(intersection);

                        if (distance < closestDistance && distance < DrawDistance)
                        {
                            nearestObject = obj;
                            closestDistance = distance;
                        }
                    }
                }

                return nearestObject;
            }

            public Vector AngleToDirectionVector(Vector baseDirection, double angle)
            {
                double rad = angle * (Math.PI / 180), cos = Math.Cos(rad), sin = Math.Sin(rad);

                Vector rot = new Vector(baseDirection.x * cos - baseDirection.z * sin, baseDirection.y, baseDirection.x * sin + baseDirection.z * cos);
                return rot;
            }

            public void SendRays(Map map, int numRays, int numVerticalRays)
            {
                double angleStep = Fov / (numRays - 1); 
                double angleStepVertical = Fov / (numVerticalRays - 1); 
                double verticalAngleStart = -Fov / 2 + angleStepVertical / 2;
                double maxAngleDeviation = Math.Tan(Fov / 2) / Math.Max(numRays, numVerticalRays); ;

                for (int j = 0; j < numVerticalRays; j++)
                {
                    double verticalRayAngle = verticalAngleStart + j * angleStepVertical;

                    for (int i = 0; i < numRays; i++)
                    {
                        double rayAngle = -Fov / 2 + i * angleStep;
                        Vector rayDirection = AngleToDirectionVector(Look_dir, rayAngle);
                        Ray ray = new Ray(position, rayDirection);

                        Point intersection = null;

                        foreach (Object obj in map.arrObj)
                        {
                            intersection = obj.Intersect(ray);
                            if (intersection != null)
                            {
                                break;
                            }
                        }

                        Vector verticalRayDirection = AngleToDirectionVector(Look_dir, verticalRayAngle);
                        Ray verticalRay = new Ray(position, verticalRayDirection);

                        Point intersectionVertical = null;

                        foreach (Object obj in map.arrObj)
                        {
                            intersectionVertical = obj.Intersect(verticalRay);
                            if (intersectionVertical != null)
                            {
                                break;
                            }
                        }

                        if (intersectionVertical != null && intersection != null)
                        {
                            // обработка пересечения
                        }
                    }
                }
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
            public Canvas(Map map, Camera camera, VectorSpace vSpace)
            {
                this.map = map;
                this.camera = camera;
                this.vSpace = vSpace;
            }
            public void Draww(Map map, Camera camera)
            {
                //метод отрисовки проекции карты map на камеру camera относительно текущего формата отрисовки.
            }
           
        }

        public class Consoles : Canvas
        {
            private ConsoleColor[] gradient = new ConsoleColor[]
            {
                ConsoleColor.Red, ConsoleColor.Green, ConsoleColor.Blue, ConsoleColor.Yellow, ConsoleColor.DarkYellow, ConsoleColor.DarkCyan

            };
            public Consoles(Map map, Camera camera, VectorSpace vSpace) : base(map, camera, vSpace) { }

            //private string dropline = "░░▒▒▓▓██ ";
            private string dropline = "@#%&*+=-:. ";



            public void Draw()
            {
                int screenWidth = Console.WindowWidth;
                int screenHeight = Console.WindowHeight;

                for (int y = 0; y < screenHeight; y++)
                {
                    for (int x = 0; x < screenWidth; x++)
                    {
                        
                        Vector dir = vSpace.Dir1 * ((double)x / screenWidth - 0.5) + vSpace.Dir2 * ((double)y / screenHeight - 0.5) + vSpace.Dir3;

                        Object obj = camera.NearestObject(map.arrObj);

                        if (obj != null)
                        {

                            Ray ray = new Ray(camera.position, dir);
                            Point intersection = obj.Intersect(ray);

                            if (intersection != null)
                            {
                                double distance = camera.DistanceTo(intersection);

                                int gradientIndex = (int)(distance / camera.DrawDistance * dropline.Length);
                                gradientIndex = Math.Min(gradientIndex, dropline.Length - 1);

                                Console.Write(dropline[gradientIndex]);
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

                    Console.WriteLine();
                }
            }


        }
        public class Angle : Point
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

            public static Point operator *(Point a, Angle angle)
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

                return new Point(newX, newY, newZ);

            }


        }

        public class Spectator : Camera
        {
            public Spectator(Point position, Vector rotation, Point Look_at, Vector Look_dir, double Fov, double DrawDistance) : base(position, rotation, Look_at, Look_dir, Fov, DrawDistance) {}

            public void MoveForward(double speed)
            {
                position += Look_at * speed;
            }

            public void MoveBackward(double speed) 
            { 
                position -= Look_at * speed;
            }


            public void Rotate(Vector direction)
            {
                rotation += direction * 2;
            }

        }

        public class Player : Camera
        {
            public Vector speed; //Скорость перемещения

            public Player(Point position, Vector rotation, Point Look_at, Vector Look_dir, double Fov, double DrawDistance) : base(position, rotation, Look_at, Look_dir, Fov, DrawDistance) 
            { 
                this.position = position;
                this.rotation = rotation;
            }

            public void Step(Vector direction) //Перемещение в указанном направлении
            {
                
            }
        }

        





    }


}

