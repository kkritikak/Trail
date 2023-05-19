///this is the change I want to see
#ifndef MPCD_SPHERE_GEOMETRY_H_
#define MPCD_SPHERE_GEOMETRY_H_

#include "BoundaryCondition.h"

#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include <cmath>

#ifdef NVCC
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#include <string>
#endif // NVCC
//comment
namespace mpcd
{
namespace detail
{
//! Sphere geometry
/*!
 * This models a fluid confined inside a sphere, centered at the origin and with radius R.
 *
 * If a particle leaves the sphere in a single simulation step, the particle is backtracked to the point on the
 * surface from which it exited the surface and then reflected according to appropriate boundary condition.
 */
class __attribute__((visibility("default"))) SphereGeometry
    {
    public:
        //! Constructor
        /*!
         * \param R confinement radius initially at 0 time
         * \param bc Boundary condition at the wall (slip or no-slip)
         * \param V velocity of spherical container (V is negative)
         * \param t is time 
         * \m_R0 is the radius of container at time t, R = R0 - V(t-t')  ,t' is initial time taking as 0 so : R = R0 - V*t
         */
        HOSTDEVICE SphereGeometry(Scalar R, Scalar V, Scalar t, boundary bc)
            : m_R(R), m_R2(R*R), m_bc(bc), m_V(V), m_V2(V*V), m_t(t), m_R0(R+V*t)
            { }

        //! Detect collision between the particle and the boundary
        /*!
         * \param pos Proposed particle position
         * \param vel Proposed particle velocity
         * \param dt Integration time remaining(time particle spent outside)
         *
         * \returns True if a collision occurred, and false otherwise
         *
         * \post The particle position \a pos is moved to the point of reflection, the velocity \a vel is updated
         *       according to the appropriate bounce back rule, and the integration time \a dt is decreased to the
         *       amount of time remaining.
         */
        HOSTDEVICE bool detectCollision(Scalar3& pos, Scalar3& vel, Scalar& dt) const
            {

            /*
             * If particle is still inside the sphere or has zero speed, no collision could have occurred and therefore
             * exit immediately.
             */
            const Scalar r2 = dot(pos,pos);
            const Scalar v2 = dot(vel,vel);
            if (r2 <= m_R0 || v2 == Scalar(0))
               {
               dt = Scalar(0);
               return false;
               }

            /*
             * Find the time remaining when the particle collided with the sphere of radius R. This time is
             * found by backtracking the position, r* = r-dt*v, and solving for dt when dot(r*,r*) = R0^2.
             * This gives a quadratic equation in dt; the smaller root is the solution.
             */

            const Scalar rv = dot(pos,vel);
            const Scalar RV = m_R0*m_V; 
            const Scalar R2 = m_R0*m_R0;
            dt = ((rv -RV) - fast::sqrt((rv-RV)*(rv-RV)-(v2-m_V2)*(r2-R2)))/(v2-m_V2);

            // backtrack the particle for time dt to get to point of contact
            pos -= vel*dt;

            // update velocity according to boundary conditions
            /*
             * Let n = r/R be the normal unit vector at the point of contact r (the particle position, which has been
             * backtracked to the surface in the previous step).
             * The perpendicular and parallel components of the velocity are:
             * v_perp = (v.n)n = (v.r/R^2)r
             * v_para = v-v_perp
             */
            if (m_bc == boundary::no_slip)
                {
                /* No-slip and no penetration requires reflection of both parallel and perpendicular components.
                 * This results in just flipping of all the velocity components.
                 */
                vel = -vel;
                }
            else if (m_bc == boundary::slip)
                {
                /*
                 * Only no-penetration condition is enforced, so only v_perp is reflected.
                 * The new velocity v' is:
                 * v' = -v_perp + v_para = v - 2*v_perp
                */
                const Scalar3 vperp = (dot(vel,pos)/m_R2)*pos;
                vel -= Scalar(2)*vperp;
                }
            return true;
            }

        //! Check if a particle is out of bounds
        /*!
         * \param pos Current particle position
         * \returns True if particle is out of bounds, and false otherwise
         */
        HOSTDEVICE bool isOutside(const Scalar3& pos) const
            {
            return dot(pos,pos) > m_R2;
            }

        //! Validate that the simulation box is large enough for the geometry
        /*!
         * \param box Global simulation box
         * \param cell_size Size of MPCD cell
         *
         * The box is large enough if the shell is padded along the radial direction, so that cells at the boundary
         * would not interact with each other via PBC.
         *
         * It would be enough to check the padding along the x,y,z directions individually as the box boundaries are
         * closest to the sphere boundary along these axes.
         */
        HOSTDEVICE bool validateBox(const BoxDim& box, Scalar cell_size) const
            {
            Scalar3 hi;
            Scalar3 lo;
            hi.x = box.getHi().x; hi.y = box.getHi().y; hi.z = box.getHi().z;
            lo.x = box.getLo().x; lo.y = box.getLo().y; lo.z = box.getLo().z;

            return ((hi.x-m_R) >= cell_size && (-lo.x-m_R) >= cell_size &&
                    (hi.y-m_R) >= cell_size && (-lo.y-m_R) >= cell_size &&
                    (hi.y-m_R) >= cell_size && (-lo.y-m_R) >= cell_size );
            }

        //! Get Sphere radius
        /*!
         * \returns confinement radius
         */
        HOSTDEVICE Scalar getR() const
            {
            return m_R0;
            }

        //! Get the wall boundary condition
        /*!
         * \returns Boundary condition at wall
         */
        HOSTDEVICE boundary getBoundaryCondition() const
            {
            return m_bc;
            }

        #ifndef NVCC
        //! Get the unique name of this geometry
        static std::string getName()
            {
            return std::string("Sphere");
            }
        #endif // NVCC

    private:
        const Scalar m_R;       //!< Sphere radius
        const Scalar m_R2;      //!< Square of sphere radius
        const boundary m_bc;    //!< Boundary condition
        const Scalar m_V;       //!<velocity of the moving spherical surface
        const Scalar m_V2;      //!<square of velocity moving surface
        const Scalar m_t;      //!<time
        const Scalar m_R0;       //!<radius of container at that time
    };

} // end namespace detail
} // end namespace mpcd

#undef HOSTDEVICE

#endif // MPCD_SPHERE_GEOMETRY_H_
