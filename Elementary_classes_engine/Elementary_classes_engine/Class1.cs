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
            public double Fov; 
            public double VFov; 
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
            public abstract double DistanceObject(Object obj);
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

            public abstract bool Intersect(Ray)
            {
                //Пересечение луча ДОБАВИТЬ
                // Заглушка для определения точки пересечения прямой с объектом
            }
            // Заглушка для определения ближайшей точки из набора
            public virtual point getNearestPoint(point[] pts)
            {
                point nearestPoint = null;
                double closestDistance = double.PositiveInfinity;
                foreach (point p in pts)
                {
                    double distance = distanceTo(point);
                    if (distance < closestDistance)
                    {
                        nearestPoint = p;
                        closestDistance = distance;
                    }
                }
                return nearestPoint;
            }
            protected double distanceTo(Point point)
            {
                double dx = point.x - position.x;
                double dy = point.y - position.y;
                double dz = point.z - position.z;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }

        }
       
        public class Plane : Object // Плоскость наследуется от объекта 
        {
            public ParametersPlane parameters;

            public Plane(point position, Vector rotation)
            {
                this.position = position;
                this.rotation = rotation;
                Parameters = new ParametersPlane();
            }

            public override point intersect(Vector V, point pt)
            {
                double denominator = parameters.CalculateValue(direction);
                if (denominator == 0) // прямая параллельна плоскости или лежит в плоскости
                {
                    if (Contains(point)) // прямая лежит в плоскости
                    {
                        return point; // любая точка прямой будет являться точкой пересечения
                    }
                    else
                    {
                        return null; // прямая параллельна плоскости и не пересекает ее
                    }
                }
                else // прямая пересекает плоскость
                {
                    double numerator = parameters.CalculateValue(position) - parameters.CalculateValue(point);
                    double distance = numerator / denominator;
                    return new Point(point.x + direction.x * distance, point.y + direction.y * distance, point.z + direction.z * distance);
                }
               
            }

            public override bool Contains(point Pt)
            {
                return parameters.CalculateValue(point); //Лежит ли точка ДОБАВИТЬ (bool) + 
            }
            
        }
        public class Parameters
        {
            //Методы преобразования коэффициентов путем поворота, перемещения и/или масштабирования
            public double aCoefficient;
            public double bCoefficient;
            public double cCoefficient;
        
            public Parameters(double a, double b, double c)
            {
                this.aCoefficient = a;
                this.bCoefficient = b;
                this.cCoefficient = c;
            }
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
            public double CalculateValue(point pt)
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
            public BoundedPlane(Vector Delta_U, Vector Delta_V, point Position, Vector Rotation) : base(Position, Rotation)
            {
                this.Delta_U = Delta_U;
                this.Delta_V = Delta_V;
            }
            public abstract bool InBoundaries(point Pt)  // Проверка координат точки на соответствие границам плоскости
            {
           
            }

        }


        public class Sphere : Object
        {
            
            private double Radius;
            private Parameters parameters;


            public Sphere(point Center, double Radius)
            {
                Position = Center;
                this.Radius = Radius;
                parameters = new Parameters;
            }
            
            public override bool Contains(point Pt)
            {
                return Math.Pow(Pt.x - Position.x, 2) + Math.Pow(Pt.y - Position.y, 2) + Math.Pow(Pt.z - Position.z, 2) == Math.Pow(Radius, 2);
                //Лежит ли точка на поверхности сферы 
            }
           
            public intersect(Vector V)
            {
                //  ДОБАВИТЬ
            }
            public static void ParametersSphere(Parameters P)
            {
                // Радиус
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
            point position;
            Vector direction;
            public Ray(point position, Vector direction) 
            {
                this.position = position;
                this.direction = direction;
            }
            public override point Intersect(Map); //метод поиска ближайшей точки пересечения с объектами на карте
            

        }
        public abstract class Canavas
        {

        }






    }


}
