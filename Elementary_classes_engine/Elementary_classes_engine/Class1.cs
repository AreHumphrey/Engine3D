using System;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Security.Cryptography.X509Certificates;
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

            public double IntersectionTo(Point point)
            {
                return Math.Sqrt((x - point.x) * (x - point.x) + (y - point.y) * (y - point.y) + (z - point.z) * (z - point.z));
            }

            public static double Distance(Point a, Point b)
            {
                return Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
            }

            public Point Subtract(Vector other) {
                return new Point(this.x - other.x, this.y - other.y, this.z - other.z );
            }

            public double GetCoordinate(int index)
            {
                if (index == 0)
                {
                    return x;
                }
                if (index == 1)
                {
                    return y;
                }
                if (index == 2) 
                { 
                    return z; 
                }
                throw new IndexOutOfRangeException("Индекс точки все диапазона");

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
            public Point pt, pt1, pt2;

            public Vector(Point pt)
            {
                this.x = pt.x;
                this.y = pt.y;
                this.z = pt.z;

            }

            public Vector(Point pt1, Point pt2)
            {
                this.pt1 = pt1;
                this.pt2 = pt2;
                x = pt2.x - pt1.x;
                y = pt2.y - pt1.y;
                z = pt2.z - pt1.z;
            }

            public Vector(double x, double y, double z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }

            public static double Lenght(Vector v1)
            {
                return Math.Sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
            }

            public Vector Normalized(Vector vector)
            {
                return vector / Lenght(vector);
            }

            public  Vector Cross(Vector a, Vector b)
            {
                double x = a.y * b.z - a.z * b.y;
                double y = a.z * b.x - a.x * b.z;
                double z = a.x * b.y - a.y * b.x;
                return new Vector(x, y, z);
            }

            public double GetCoordinate(int index)
            {
                if (index == 0) return x;
                if (index == 1) return y;
                if (index == 2) return z;
                throw new IndexOutOfRangeException("Индекс точки все диапазона");
            }

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
                double vectorsLen = Vector.Lenght(vect);
                return vect / vectorsLen;
            }

        }



        public class Camera
        {
            public double width = 70;
            public double height = 50;
            public double fov;
            public double vFov; 
            public Vector lookAt;
            public Vector lookDir;
            public Point position;
            public Vector rotation;
            public double drawDistance;

            public Camera(Point position, Vector rotation, Vector lookAt, Vector lookDir, double fov, double drawDistance)
            {
                this.rotation = rotation;
                this.lookAt = lookAt;
                this.lookDir = lookDir;
                this.fov = (fov / 180 * Math.PI) / 2;
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

            public Object NearestObject(Object[] objects, Vector rayDirection)
            {
                Object nearestObject = null; double closestDistance = double.PositiveInfinity;

                foreach (Object obj in objects)
                {
                    Ray ray = new Ray(position, rayDirection);
                    Point intersection = obj.Intersect(ray);
                    if (intersection != null)
                    {
                        double distance = DistanceTo(intersection);
                        if (distance < closestDistance && distance < drawDistance)
                        {
                            bool isCloser = true;
                            foreach (Object other in objects)
                            {
                                if (other != obj)
                                {
                                    Ray otherRay = new Ray(position, rayDirection);
                                    Point otherIntersection = other.Intersect(otherRay);
                                    if (otherIntersection != null)
                                    {
                                        double otherDistance = DistanceTo(otherIntersection);
                                        if (otherDistance < distance)
                                        {
                                            isCloser = false;
                                            break;
                                        }
                                    }
                                }
                            }

                            if (isCloser)
                            {
                                nearestObject = obj;
                                closestDistance = distance;
                            }
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

            public List<Ray> SendRays(Map map)
            {
                List<Ray> rays = new List<Ray>();
                double rayAngle = fov / width, vRayAngle = vFov / height; 

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
            public Vector rotation;

            public abstract bool Contains(Point p);

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
                double dx = point.x - position.x, dy = point.y - position.y, dz = point.z - position.z;

                return Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }

            public Point GetCollisionPoint(Point playerPosition, Vector playerDirection)
            {
                Ray playerRay = new Ray(playerPosition, playerDirection);
                Point intersectionPoint = Intersect(playerRay);

                return intersectionPoint;
            }

        }



        public class Plane : Object
        {
            public Parameters parameters;

            public Plane(Point position, Vector rotation)
            {
                this.position = position;
                this.rotation = rotation;
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



        public class BoundedPlane : Plane
        {
            private double dx, dy, dz; 

            public BoundedPlane(Point position, Vector rotation, double dx, double dy, double dz) : base(position, rotation)
            {
                this.dx = dx;
                this.dy = dy;
                this.dz = dz;
            }

            public override bool Contains(Point point)
            {
                bool contains = base.Contains(point);

                if (contains)
                {
                    double xDiff = point.x - position.x, yDiff = point.y - position.y, zDiff = point.z - position.z;
                    contains = (Math.Abs(xDiff) <= dx / 2) && (Math.Abs(yDiff) <= dy / 2) && (Math.Abs(zDiff) <= dz / 2);
                }

                return contains;
            }

            public override Point Intersect(Ray ray)
            {
                return base.Intersect(ray); 
            }

            public override Point NearestPoint(Point[] points)
            {
                Point[] boundedPoints = points.Where(p => Contains(p)).ToArray(); 

                return base.NearestPoint(boundedPoints); 
            }
        }



        public class Parameters
        {
            public double aCoefficient, bCoefficient, cCoefficient;

            public void Rotate(double angle) 
            {
                double cosA = Math.Cos(angle);
                double sinA = Math.Sin(angle);

                double newACoefficient = aCoefficient * cosA - bCoefficient * sinA;
                double newBCoefficient = bCoefficient * cosA + aCoefficient * sinA;

                aCoefficient = newACoefficient;
                bCoefficient = newBCoefficient;
            }

            public void Move(double x, double y, double z) 
            {
                aCoefficient += x;
                bCoefficient += y;
                cCoefficient += z;
            }

            public void Scale(double xFactor, double yFactor, double zFactor)  
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
            public double c = 1;
            public double d;  
            
            public ParametersPlane()
            {
                a = 0;
                b = 0;
                c = 1;
                d = 0;
            }
            
        }



        public class Parallelepiped : Object
        {
            private ParametersParallelepiped parameters;

            public Parallelepiped(Point center, ParametersParallelepiped parameters)
            {
                position = center;
                this.parameters = parameters;
            }

            public override bool Contains(Point point)
            {
                double distanceX = Math.Abs(point.x - position.x);
                double distanceY = Math.Abs(point.y - position.y);
                double distanceZ = Math.Abs(point.z - position.z);

                if (distanceX > parameters.width / 2 || distanceY > parameters.height / 2 || distanceZ > parameters.depth / 2)
                {
                    return false;
                }

                return true;
            }

            public override Point Intersect(Ray ray)
            {
                double tmin = double.MinValue, tmax = double.MaxValue;
                Point p = position;
                double[] bounds = { p.x - parameters.width / 2, p.x + parameters.width / 2, p.y - parameters.height / 2, p.y + parameters.height / 2, p.z - parameters.depth / 2, p.z + parameters.depth / 2 }; 

                for (int i = 0; i < bounds.Length; i += 2)
                {
                    double t0 = (bounds[i] - ray.position.GetCoordinate(i / 2)) / ray.direction.GetCoordinate(i / 2); 
                    double t1 = (bounds[i + 1] - ray.position.GetCoordinate(i / 2)) / ray.direction.GetCoordinate(i / 2);

                    if (t0 > t1)
                    {
                        double temp = t0;
                        t0 = t1;
                        t1 = temp;
                    }

                    if (t0 > tmin)
                    {
                        tmin = t0;
                    }

                    if (t1 < tmax)
                    {
                        tmax = t1;
                    }

                    if (tmin > tmax)
                    {
                        return null;
                    }

                    if (tmax < 0)
                    {
                        return null;
                    }
                }

                double t = tmin > 0 ? tmin : tmax;

                Point point = new Point(ray.position.x + t * ray.direction.x, ray.position.y + t * ray.direction.y, ray.position.z + t * ray.direction.z);

                if (Contains(point))
                {
                    return point;
                }

                return null;
            }
        }



        public class ParametersParallelepiped : Parameters
        {
            public double width;
            public double height;
            public double depth;

            public ParametersParallelepiped(double x, double y, double z, double width, double height, double depth)
            {
                this.width = width;
                this.height = height;
                this.depth = depth;
                this.aCoefficient = x;
                this.bCoefficient = y;
                this.cCoefficient = z;
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
                    Point point = obj.Intersect(this); 
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
            public void Draw() { }
          
        }




        public class Consoles : Canvas
        {
            public Consoles(Map map, Camera camera, VectorSpace vSpace) : base(map, camera, vSpace) { }

            private string dropline = "@#*+=-^:. ";

            public void Draw()
            {
                int screenWidth = Console.WindowWidth - 45;
                int screenHeight = Console.WindowHeight;

                for (int y = 0; y < screenHeight; y++)
                {
                    for (int x = 0; x < screenWidth; x++)
                    {
                        Vector dir = vSpace.dir1 * ((double)x / screenWidth - 0.5) + vSpace.dir2 * ((double)y / screenHeight - 0.5) + vSpace.dir3;
                        Object nearestObject = camera.NearestObject(map.arrObj, dir);

                        if (nearestObject != null)
                        {
                            Ray ray = new Ray(camera.position, dir);
                            Point intersection = nearestObject.Intersect(ray);
                            if (intersection != null)
                            {
                                double distance = camera.DistanceTo(intersection);
                                double gradientIndex = (distance / 100.0) * dropline.Length;
                                for (int i = 1; i < nearestObject.DistanceTo(camera.position); i++)
                                {
                                    gradientIndex += 0.05; 
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
                double newYZ = a.yz + b.yz, newXZ = a.xz + b.xz, newXY = a.xy + b.xy;

                return new Angle(0, 0, 0, newYZ, newXZ, newXY);
            }

            public static Angle operator -(Angle a)
            {
                double newYz = -a.yz, newXz = -a.xz, newXy = -a.xy;

                return new Angle(0, 0, 0, newYz, newXz, newXy);
            }

            public static Vector operator *(Point a, Angle angle)
            {
                double yzRad = angle.yz * Math.PI / 180, xzRad = angle.xz * Math.PI / 180,  xyRad = angle.xy * Math.PI / 180;

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
            public Map map;

            public Player( Camera camera, double drawDistance, double collisionRadius, Map map)
            {
                this.camera = camera;
                this.map = map;
            }

            public Camera MoveForward(double distance)
            {
                double dx = distance * camera.lookAt.x;
                double dz = distance * camera.lookAt.z;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);
                Point collisionPt = CheckCollision(newPosition, playerDirection);

                if(collisionPt != newPosition)
                { 
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }
                return camera;
            }

            public Camera MoveBack(double distance)
            {
                double dx = distance * camera.lookAt.x, dz = distance * camera.lookAt.z;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x - dx, camera.position.y, camera.position.z - dz);
                Point collisionPt = CheckCollision(newPosition, playerDirection);
                if (collisionPt != newPosition)
                {
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }
                return camera;

            }

            public Camera MoveUp(double distance)
            {
                double dx = distance * camera.lookDir.x, dy = distance * camera.lookDir.y, dz = distance * camera.lookDir.z;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x, camera.position.y - distance, camera.position.z);
                Point collisionPt = CheckCollision(newPosition, playerDirection);

                if (collisionPt != newPosition)
                {
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }

                return camera;
                
            }

            public Camera MoveDown(double distance)
            {
                double dx = distance * camera.lookDir.x, dy = distance * camera.lookDir.y, dz = distance * camera.lookDir.z;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x, camera.position.y +  distance, camera.position.z);
                Point collisionPt = CheckCollision(newPosition, playerDirection);

                if (collisionPt != newPosition)
                {
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }

                return camera;

            }

            public Camera StrafeLeft(double distance)
            {
                double dx = -distance * camera.lookAt.z, dz = -distance * camera.lookAt.x;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);
                Point collisionPt = CheckCollision(newPosition, playerDirection);

                if (collisionPt != newPosition)
                {
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }

                return camera;
            }

            public Camera StrafeRight(double distance)
            {
                double dx = -distance * camera.lookAt.z, dz = distance * camera.lookAt.x;
                Vector playerDirection = camera.lookAt;
                Point newPosition = new Point(camera.position.x - dx, camera.position.y, camera.position.z - dz);
                Point collisionPt = CheckCollision(newPosition, playerDirection);
                if ((collisionPt.x < newPosition.x && collisionPt.y < newPosition.y && collisionPt.z < newPosition.z))
                {
                    camera.position = collisionPt;
                }
                else
                {
                    camera.position = newPosition;
                }

                return camera;
            }

            public Camera Rotate(double degX, double degY, double degZ)
            {
                degX = degX % 360; degY = degY % 360; degZ = degZ % 360;

                camera.position = new Point(camera.position.x + degX, camera.position.y + degY, camera.position.z + degZ);
                double radX = camera.position.x * Math.PI / 180, radY = camera.position.y * Math.PI / 180, radZ = camera.position.z * Math.PI / 180;
                camera.lookDir = new Vector(Math.Sin(radY) * Math.Cos(radX), Math.Sin(radX), Math.Cos(radY) * Math.Cos(radX));
                camera.position = new Point(camera.position.x, camera.position.y, camera.position.z);
                camera.lookAt.x = camera.position.x + camera.lookDir.x;
                camera.lookAt.y = camera.position.y + camera.lookDir.y;
                camera.lookAt.z = camera.position.z + camera.lookDir.z;

                return camera;
            }

            public Point CheckCollision(Point newPosition, Vector playerDirection)
            {
                foreach (Object obj in map.arrObj)
                {
                    if (obj.Contains(newPosition))
                    {
                        Point collisionPoint = obj.GetCollisionPoint(newPosition, playerDirection);
                        return collisionPoint;
                    }
                }

                return newPosition;
            }


        }




        public class Spectator
        {
            public Camera camera;
            public Point position;
            public Vector lookAt;
            public double dx;

            public Spectator(Camera camera, double drawDistance)
            {
                this.camera = camera;
            }

            public Camera MoveBack(double distance) 
            {
                double dx = distance * camera.lookDir.x, dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x - dx, camera.position.y, camera.position.z - dz);

                return camera;
            }

            public Camera MoveForward(double distance)
            {
                double dx = distance * camera.lookDir.x, dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);

                return camera;
            }

            public Camera MoveUp(double distance)
            {
                double dx = distance * camera.lookDir.x, dy = distance * camera.lookDir.y, dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x + dx, camera.position.y + dy, camera.position.z);

                return camera;
            }
            public Camera MoveDown(double distance)
            {
                double dx = distance * camera.lookDir.x, dy = distance * camera.lookDir.y, dz = distance * camera.lookDir.z;
                camera.position = new Point(camera.position.x, camera.position.y + distance, camera.position.z);
                return camera;
               
            }


            public Camera StrafeLeft(double distance)
            {
                double dx = distance * camera.lookDir.z, dz = -distance * camera.lookDir.x;
                camera.position = new Point(camera.position.x + dx, camera.position.y, camera.position.z + dz);

                return camera;
            }

            public Camera StrafeRight(double distance)
            {
                double dx = -distance * camera.lookDir.z, dz = distance * camera.lookDir.x;
                camera.position = new Point(camera.position.x - dx, camera.position.y, camera.position.z - dz);

                return camera;
            }

            public Camera Rotate(double degX, double degY, double degZ)
            {
                camera.position = new Point(camera.position.x + degX, camera.position.y + degY, camera.position.z + degZ);
                double radX = camera.position.x * Math.PI / 180, radY = camera.position.y * Math.PI / 180, radZ = camera.position.z * Math.PI / 180;
                camera.lookDir = new Vector(Math.Sin(radY) * Math.Cos(radX), Math.Sin(radX), Math.Cos(radY) * Math.Cos(radX));
                Vector up = new Vector(Math.Sin(radY + Math.PI / 2) * Math.Cos(radX), Math.Sin(radX + Math.PI / 2), Math.Cos(radY + Math.PI / 2) * Math.Cos(radX));
                Vector right = up.Cross(camera.lookDir, up);
                camera.position = new Point(camera.position.x, camera.position.y, camera.position.z);
                camera.lookAt.x = camera.position.x + camera.lookDir.x;
                camera.lookAt.y = camera.position.y + camera.lookDir.y;
                camera.lookAt.z = camera.position.z + camera.lookDir.z;

                return camera;
            }
           
        }



        public class Event
        {
            Dictionary<string, Action<object[]>> events = new Dictionary<string, Action<object[]>>();

            public void Append(string name)
            {
                events.Add(name, null);
            }

            public void Handle(string name, Action<object[]> callfunction)
            {
                events[name] = callfunction;
            }

            public void Trigger(string name, params object[] args)
            {
                if (events.ContainsKey(name))
                {
                    events[name]?.Invoke(args);
                }
            }
        }



        public class Trigger
        {
            Event eventSystem;

            public Trigger(Event eventSystem)
            {
                this.eventSystem = eventSystem;
            }

            public void Triggers(string name, params object[] args)
            {
                eventSystem.Trigger(name, args);
            }
        }
        
    }


}

