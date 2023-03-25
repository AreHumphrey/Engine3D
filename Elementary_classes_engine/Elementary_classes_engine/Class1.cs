using System;
using System.ComponentModel;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Threading;

namespace Elementary_classes_engine
{
    public class Program
    {
        public class point
        {
            
            public double x { get; private set; }
            public double y { get; private set; }
            public double z { get; private set; }
            public point(double x, double y, double z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }
            //Перегрузка операторов ( + - * / ) + расстояние
            public static double distance(point a, point b)
            {
                return Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
            }
            public static point operator +(point a, point b)
            {
                return new point(a.x + b.x, a.y + b.y, a.z + b.z);
            }
            public static point operator -(point a, point b)
            {
                return new point(a.x - b.x, a.y - b.y, a.z - b.z);
            }
            public static point operator *(point a, point b)
            {
                return new point(a.x * b.x, a.y * b.y, a.z * b.z);
            }
            public static point operator /(point a, point b)
            {
                return new point(a.x / b.x, a.y / b.y, a.z / b.z);
            }
           

            public static point operator *(point a, double Element)
            {
                return new point(a.x * Element, a.y * Element, a.z * Element);
            }
            public static point operator *(double Element, point b)
            {
                return new point(b.x * Element, b.y * Element, b.z * Element);
            }
            public static point operator /(point a, double Element)
            {
                return new point(Element / a.x, Element / a.y, Element / a.z);
            }
        }
        public class Vector
        {
            public point pt;

            public Vector(point pt)
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
                return new Vector(new point(a.pt.x + b.pt.x, a.pt.y + b.pt.y, a.pt.z + b.pt.z));
            }
            public static Vector operator -(Vector a, Vector b)
            {
                return new Vector(new point(a.pt.x - b.pt.x, a.pt.y - b.pt.y, a.pt.z - b.pt.z));
            }
            public static Vector operator *(Vector a, Vector b)
            {
                return new Vector(new point(a.pt.x * b.pt.x, a.pt.y * b.pt.y, a.pt.z * b.pt.z));
            }
            public static Vector operator ^(Vector a, Vector b)
            {
                return new Vector(new point(a.pt.y * b.pt.z - a.pt.z * b.pt.z, a.pt.z * b.pt.x - a.pt.x * b.pt.z, a.pt.x * b.pt.y - a.pt.y * b.pt.x));
            }
            public static Vector operator /(Vector a, Vector b)
            {
                return new Vector(new point(a.pt.x / b.pt.x, a.pt.y / b.pt.y, a.pt.z / b.pt.z));
            }
            public static Vector operator *(double Element, Vector b)
            {
                return new Vector(new point(Element * b.pt.x, Element * b.pt.y, Element * b.pt.z));
            }
            public static Vector operator *(Vector a, double Element)
            {
                return new Vector(new point(Element * a.pt.x, Element * a.pt.y, Element * a.pt.z));
            }
            public static Vector operator /(Vector a, double Element)
            {
                return new Vector(new point(a.pt.x / Element, a.pt.y / Element, a.pt.z / Element));
            }
            public static Vector operator /(double Element, Vector b)
            {
                return new Vector(new point(Element / b.pt.x, Element / b.pt.y, Element / b.pt.z));
            }


        }
        public class VectorSpace
        {
            public point Initpt;
            public Vector Dir1;
            public Vector Dir2;
            public Vector Dir3;

            public VectorSpace(point InitPt, Vector Dir1, Vector Dir2, Vector Dir3)
            {
                this.Initpt = InitPt;
                this.Dir1 = VectorBasis(Dir1);
                this.Dir2 = VectorBasis(Dir1);
                this.Dir3 = VectorBasis(Dir1);
            }

            public static Vector VectorBasis(Vector vect)
            {
                double VectorsLen = Vector.LenVector(vect);
                return vect / VectorsLen;
            }
        }
        public abstract class Camera
        {
            public double Width = 70;
            public double Height = 50;
            public double Fov; //d +
            public double VFov; //d + 
            public point Look_at;
            public double Look_dir;
            public double DrawDistance;
            public Camera(point Pos, point Look_at, double Look_dir, double Fov, double DrawDistance, double vFov)
            {
                this.Look_at = Look_at;
                this.Look_dir = Look_dir;
                this.Fov = Fov;
                this.VFov = vFov;
                vFov = Fov * (Height / Width);
                this.DrawDistance = DrawDistance;

            }
            public abstract double Func();// Пустить во все стороны от камеры лучи (необходим подсчет)

        }
        public static double Sum(double a, double b)
        {
            return a + b;
        }
        public abstract class Object
        {
            public point Position;
            public Vector Rotation; //пока оставить
            public Object(point Position, Vector Rotation)
            {
                this.Position = Position;
                this.Rotation = Rotation;
            }

            public abstract bool Contains(point Pt); //Лежит ли точка ДОБАВИТЬ (bool)
            
            public static void Intersect(Vector V)
            {
                 //Пересечение луча ДОБАВИТЬ
            }

          
        }
       
        public class Plane : Object // Плоскость наследуется от объекта 
        {
            public Plane(point Position, Vector Rotation) : base(Position, Rotation){}
            public static void contains(point Pt)
            {
                //Лежит ли точка ДОБАВИТЬ (bool)
            }
            public static void intersect(Vector V)
            {
                // Точка Пересечения камеры и объекта ДОБАВИТЬ
            }

            public override bool Contains(point Pt)
            {
                throw new NotImplementedException();
            }
        }
        public class Parameters
        {
            //методы преобразования коэффициентов путем поворота, перемещения и/или масштабирования
            public double aCoefficient;
            public double bCoefficient;
            public double cCoefficient;
            public Parameters(double a, double b, double c)
            {
                this.aCoefficient = a;
                this.bCoefficient = b;
                this.cCoefficient = c;
            }
        }
        public class BoundedPlane : Plane //  Главная плоскость
        {
            Vector Delta_U;
            Vector Delta_V;
            public BoundedPlane(Vector Delta_U, Vector Delta_V, point Position, Vector Rotation) : base(Position, Rotation)
            {
                this.Delta_U = Delta_U;
                this.Delta_V = Delta_V;
            }
            public static bool InBoundaries(point Pt)  // Проверка координат точки на соответствие границам плоскости
            {
                return true;
            }

        }


        public class Sphere : Object
        {
            public double Radius;
            
          

            public Sphere(point Center, double Radius, point Position, Vector Rotation) : base(Position, Rotation)
            {
                Position = Center;
                this.Radius = Radius;
               
            }
            
            public override bool Contains(point Pt)
            {
                return Math.Pow(Pt.x - Position.x, 2) + Math.Pow(Pt.y - Position.y, 2) + Math.Pow(Pt.z - Position.z, 2) == Math.Pow(Radius, 2);
                //Лежит ли точка 
            }
           
            public static void intersect(Vector V)
            {
                //  ДОБАВИТЬ
            }
            public static void ParametersSphere(Parameters P)
            {
                // Радиус
            }
        }


         
    }


}
