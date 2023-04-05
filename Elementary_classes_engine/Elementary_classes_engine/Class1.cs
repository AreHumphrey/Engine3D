using NPOI.SS.Formula.Functions;
using System;
using System.ComponentModel;
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
            public static Point operator /(Point a, double Element)
            {
                return new Point(Element / a.x, Element / a.y, Element / a.z);
            }
        }
        public class Vector
        {
            public Point pt;

            public Vector(Point pt)
            {
                this.pt = pt;

            }

            public static double LenVector(Vector v1)
            {
                return Math.Sqrt(v1.pt.x * v1.pt.x + v1.pt.y * v1.pt.y + v1.pt.z * v1.pt.z);
            }


            //Перегрузка операторов ( + - * / ) + расстояние
            public static Vector operator +(Vector a, Vector b)

            {
                return new Vector(new Point(a.pt.x + b.pt.x, a.pt.y + b.pt.y, a.pt.z + b.pt.z));
            }
            public static Vector operator -(Vector a, Vector b)
            {
                return new Vector(new Point(a.pt.x - b.pt.x, a.pt.y - b.pt.y, a.pt.z - b.pt.z));
            }
            public static Vector operator *(Vector a, Vector b)
            {
                return new Vector(new Point(a.pt.x * b.pt.x, a.pt.y * b.pt.y, a.pt.z * b.pt.z));
            }
            public static Vector operator ^(Vector a, Vector b)
            {
                return new Vector(new Point(a.pt.y * b.pt.z - a.pt.z * b.pt.z, a.pt.z * b.pt.x - a.pt.x * b.pt.z, a.pt.x * b.pt.y - a.pt.y * b.pt.x));
            }
            public static Vector operator /(Vector a, Vector b)
            {
                return new Vector(new Point(a.pt.x / b.pt.x, a.pt.y / b.pt.y, a.pt.z / b.pt.z));
            }
            public static Vector operator *(double Element, Vector b)
            {
                return new Vector(new Point(Element * b.pt.x, Element * b.pt.y, Element * b.pt.z));
            }
            public static Vector operator *(Vector a, double Element)
            {
                return new Vector(new Point(Element * a.pt.x, Element * a.pt.y, Element * a.pt.z));
            }
            public static Vector operator /(Vector a, double element)
            {
                return new Vector(new Point(a.pt.x / element, a.pt.y / element, a.pt.z / element));
            }
            public static Vector operator /(double Element, Vector b)
            {
                return new Vector(new Point(Element / b.pt.x, Element / b.pt.y, Element / b.pt.z));
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
                this.Dir2 = VectorBasis(Dir1);
                this.Dir3 = VectorBasis(Dir1);
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
            public double VFov; 
            public Point Look_at;
            public double Look_dir;
            public double DrawDistance;
            public Camera(Point Pos, Point Look_at, double Look_dir, double Fov, double DrawDistance, double vFov)
            {
                this.Look_at = Look_at;
                this.Look_dir = Look_dir;
                this.Fov = Fov;
                this.VFov = vFov;
                vFov = Fov * (Height / Width);
                this.DrawDistance = DrawDistance;

            }
            public static void Func(Vector v)
            {
                // Пустить во все стороны от камеры лучи (необходим подсчет)
            }
            public void sendRays(Map map) { }
        }
        public static double Sum(double a, double b)
        {
            return a + b;
        }
        public class Object
        {
            public Point Position;
            public Vector Rotation; //пока оставить
            public void Contains(Point Pt)
            {
                //Лежит ли точка ДОБАВИТЬ (bool)
            }

            public  void Intersect(Ray ray)
            {

            }

            //Пересечение луча ДОБАВИТЬ
            // Заглушка для определения точки пересечения прямой с объектом

            // Заглушка для определения ближайшей точки из набора
            public  void getNearestPoint(Point[] pts)
            {

            }
            
            public void distanceTo(Point point)
            {

            }
            
        }
       
        public class Plane : Object // Плоскость наследуется от объекта 
        {
            public ParametersPlane parameters;

            public Plane(Point position, Vector rotation)
            {
                this.Position = position;
                this.Rotation = rotation;
                parameters = new ParametersPlane();
            }

            public  void Intersect(Vector V, Point pt)
            {
               
            }

            public double Contains(Point pt)
            {
                return parameters.CalculateValue(pt); //Лежит ли точка ДОБАВИТЬ (bool) + 
            }
            
        }
        public class Parameters
        {
            //Методы преобразования коэффициентов путем поворота, перемещения и/или масштабирования
            public double aCoefficient;
            public double bCoefficient;
            public double cCoefficient;
        
            public void Rotate(double angle)
            {
                // Преобразовать коэффициенты уравнения при помощи поворота
            }

            public void Move(double x, double y, double z)
            {
                // Преобразовать коэффициенты уравнения при помощи перемещения
            }

            public void Scale(double xFactor, double yFactor, double zFactor)
            {
                // Преобразовать коэффициенты уравнения при помощи масштабирования
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
            public double CalculateValue(Point pt)
            {
                return a * pt.x + b * pt.y + c * pt.z + d;
            }
            public double CalculateValue(Vector v)
            {
                return a * v.pt.x + b * v.pt.y + c * v.pt.z;
            }
        }
        public class BoundedPlane : Plane //  Главная плоскость
        {
            Vector Delta_U;
            Vector Delta_V;
            public BoundedPlane(Vector Delta_U, Vector Delta_V, Point Position, Vector Rotation) : base(Position, Rotation)
            {
                this.Delta_U = Delta_U;
                this.Delta_V = Delta_V;
            }
            public void InBoundaries(Point Pt)  // Проверка координат точки на соответствие границам плоскости
            {
           
            }

        }


        public class Sphere : Object
        {
            
            private double Radius;
            private Parameters parameters;


            public Sphere(Point Center, double Radius)
            {
                Position = Center;
                this.Radius = Radius;
                parameters = new Parameters();
            }
            
            public bool Contains(Point Pt)
            {
                return Math.Pow(Pt.x - Position.x, 2) + Math.Pow(Pt.y - Position.y, 2) + Math.Pow(Pt.z - Position.z, 2) == Math.Pow(Radius, 2);
                //Лежит ли точка на поверхности сферы 
            }
           
            public void intersect(Vector V)
            {
                //  ДОБАВИТЬ
            }
            


        }
        public class ParametersSphere : Sphere
        {
            public ParametersSphere(Point Center, double Radius) : base(Center, Radius)
            {
            }
        }
        public class Map
        {
            private List<T> arrObjects;
            public Map()
            {
                arrObjects = new List<T>();
            }
            public void Append(T obj)
            {
                arrObjects.Add(obj);
            }
            public void Remove(T obj)
            {
                arrObjects.Remove(obj);
            }

        }
        public class Ray
        {
            Point position;
            Vector direction;
            public Ray(Point position, Vector direction) 
            {
                this.position = position;
                this.direction = direction;
            }
            public void NearestObject(Map map) 
            {
                //метод поиска ближайшей точки пересечения с объектами на карте
            }


        }
        public class Canvas
        {
            public Map map;
            public Camera camera;
            public VectorSpace vspace;
            public void Draw(Map map, Camera camera)
            {
                //метод отрисовки проекции карты map на камеру camera относительно текущего формата отрисовки.
            }

        }

        public abstract class Console : Canvas
        {
            private ConsoleColor[] gradient = new ConsoleColor[]
            {
                ConsoleColor.Red, ConsoleColor.Green, ConsoleColor.Blue, ConsoleColor.Yellow, ConsoleColor.DarkYellow, ConsoleColor.DarkCyan

            };
           
            public void Draw(Map map, Camera camera)
            {
                //Требуется перегрузить метод draw для отрисовки конкретно в консоль;

            }
            double[] graphics = {};
            public void DrawDistant(double[] graphics)
            {
                //в зависимости от дальности объекта на экране.
            }

            char[] dropline = { ' ', '#', '%', ':', ';', ',', '.' };
        }







    }


}
