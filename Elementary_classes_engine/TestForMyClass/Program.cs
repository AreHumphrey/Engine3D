
using static Elementary_classes_engine.Program;

namespace Tests_for_Class
{
    public class Program
    {
        static void Main()
        {
            point Pt1 = new point(1, 2, 5);
            point Pt2 = new point(-1, 4, 6);
            point Pt = new point(0, 0, -1);
            Console.WriteLine($"Pt1: ({Pt1.x}, {Pt1.y}, {Pt1.z}) \nPt2: ({Pt2.x}, {Pt2.y}, {Pt2.z})");

            point Pt0 = Pt + Pt2; // (0, 0, -1) + (-1, 4, 6) Pt Pt2
            Console.WriteLine($"Pt0: ({Pt0.x}, {Pt0.y}, {Pt0.z})"); // (-1, 4, 5) Pt0
            point Pt3 = Pt0 - Pt1; // (-1, 4, 5) - (1, 2, 5)
            Console.WriteLine($"Pt3: ({Pt3.x}, {Pt3.y}, {Pt3.z})"); //(-2, 2, 0) Pt3
            double d = 9;
            point Pt4 = Pt3 * d;
            Console.WriteLine($"Pt4: ({Pt4.x}, {Pt4.y}, {Pt4.z})"); //Pt4 (-18, 18, 0)
            point Pt5 = d * Pt3;
            Console.WriteLine($"Pt5: ({Pt5.x}, {Pt5.y}, {Pt5.z})"); //Pt5 (-18, 18, 0)
            point Pt6 = new point(2, 3, -2);
            point Pt7 = new point(5, 5, -25);
            Vector V1 = new Vector(Pt6);
            Vector V2 = new Vector(Pt7);
            Vector V3 = V1 + V2; //(2, 3, -2) + (5, 5, -25)
            Vector V4 = V1 - V2; //(2, 3, -2) - (5, 5, -25)
            Vector V5 = V1 * d; //(2, 3, -2) * 9
            Console.WriteLine($"V3: ({V3.pt.x}, {V3.pt.y}, {V3.pt.z})"); //V3: (7, 8, -27)
            Console.WriteLine($"V4: ({V4.pt.x}, {V4.pt.y}, {V4.pt.z})"); //V4: (-3, -2, 23)
            Console.WriteLine($"V5: ({V5.pt.x}, {V5.pt.y}, {V5.pt.z})"); //V5: (18, 27, -18)
            Vector V6 = V1 ^ V2;
            Console.WriteLine($"V6: ({V6.pt.x}, {V6.pt.y}, {V6.pt.z})"); //V6: (-125, 40, -5)
            double LenOfVectorV6 = Vector.LenVector(V6); //(Значение с онлайн калькулятора ≈ 131.33925536563697) мое Len: 131,33925536563697
            double LenOfVectorV5= Vector.LenVector(V5); // Len: 37,107950630558946 (≈ 37.107950630558946)
            double LenOfVectorV4 = Vector.LenVector(V4); //Len: 23,280893453645632 (≈ 23.280893453645632)
            double LenOfVectorV3 = Vector.LenVector(V3); //Len: 29,017236257093817 (≈ 29.017236257093817)
            Console.WriteLine($"Len V6: {LenOfVectorV6}");
            Console.WriteLine($"Len V5: {LenOfVectorV5}");
            Console.WriteLine($"Len V4: {LenOfVectorV4}");
            Console.WriteLine($"Len V3: {LenOfVectorV3}");




        }
    }
}

