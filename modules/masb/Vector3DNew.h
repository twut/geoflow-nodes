#include"amp.h"
#include"amp_math.h"
#include "Vector3D.h"
#pragma once

class Vector3DNew
{

public:
    float x;
    float y;
    float z;

    static float PointToPointDis(Vector3DNew p1, Vector3DNew p2) restrict(cpu, amp) {
        float dis = concurrency::precise_math::sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
        return dis;
    }

    static int inline GetIntersection(float fDst1, float fDst2, Vector3DNew P1, Vector3DNew P2, Vector3DNew &Hit)restrict(amp, cpu) {
        if ((fDst1 * fDst2) >= 0.0f) return 0;
        if (fDst1 == fDst2) return 0;
        Hit.x = P1.x + (P2.x - P1.x) * (-fDst1 / (fDst2 - fDst1));
        Hit.y = P1.y + (P2.y - P1.y) * (-fDst1 / (fDst2 - fDst1));
        Hit.z = P1.z + (P2.z - P1.z) * (-fDst1 / (fDst2 - fDst1));
        return 1;
    }

    static int inline InBox(Vector3DNew Hit, Vector3DNew B1, Vector3DNew B2, const int Axis)restrict(amp, cpu) {
        if (Axis == 1 && Hit.z > B1.z && Hit.z < B2.z && Hit.y > B1.y && Hit.y < B2.y) return 1;
        if (Axis == 2 && Hit.z > B1.z && Hit.z < B2.z && Hit.x > B1.x && Hit.x < B2.x) return 1;
        if (Axis == 3 && Hit.x > B1.x && Hit.x < B2.x && Hit.y > B1.y && Hit.y < B2.y) return 1;
        return 0;
    }
    static int CheckLineBox(Vector3DNew B1, Vector3DNew B2, Vector3DNew L1, Vector3DNew L2, Vector3DNew &Hit)restrict(amp, cpu)
    {
        if (L2.x < B1.x && L1.x < B1.x) return false;
        if (L2.x > B2.x && L1.x > B2.x) return false;
        if (L2.y < B1.y && L1.y < B1.y) return false;
        if (L2.y > B2.y && L1.y > B2.y) return false;
        if (L2.z < B1.z && L1.z < B1.z) return false;
        if (L2.z > B2.z && L1.z > B2.z) return false;
        if (L1.x > B1.x && L1.x < B2.x &&
            L1.y > B1.y && L1.y < B2.y &&
            L1.z > B1.z && L1.z < B2.z)
        {
            Hit = L1;
            return true;
        }
        if ((GetIntersection(L1.x - B1.x, L2.x - B1.x, L1, L2, Hit) && InBox(Hit, B1, B2, 1))
            || (GetIntersection(L1.y - B1.y, L2.y - B1.y, L1, L2, Hit) && InBox(Hit, B1, B2, 2))
            || (GetIntersection(L1.z - B1.z, L2.z - B1.z, L1, L2, Hit) && InBox(Hit, B1, B2, 3))
            || (GetIntersection(L1.x - B2.x, L2.x - B2.x, L1, L2, Hit) && InBox(Hit, B1, B2, 1))
            || (GetIntersection(L1.y - B2.y, L2.y - B2.y, L1, L2, Hit) && InBox(Hit, B1, B2, 2))
            || (GetIntersection(L1.z - B2.z, L2.z - B2.z, L1, L2, Hit) && InBox(Hit, B1, B2, 3)))
            return true;

        return false;
    }

    Vector3DNew() restrict(amp, cpu)
    {
    }
    Vector3DNew(float xx, float yy, float zz) restrict(amp, cpu)
    {
        x = xx;
        y = yy;
        z = zz;

    }
    Vector3DNew(Vector3D v1) restrict(amp, cpu)
    {
        x = v1.x;
        y = v1.y;
        z = v1.z;
    }
    Vector3DNew &operator = (const Vector3D &v1) restrict(amp, cpu)
    {
        x = v1.x;
        y = v1.y;
        z = v1.z;
    }
};


