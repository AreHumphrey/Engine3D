using System;
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
        public class Camera
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
        }
        public class Object
        {
            public point Position;
            public Vector Rotation; //пока оставить
            public Object(point Position, Vector Rotation)
            {
                this.Position = Position;
                this.Rotation = Rotation;
            }


        }



    }


}