#ifndef ADLER_PERLIN
#define ADLER_PERLIN

#include <cmath>
#include <boost/random/variate_generator.hpp>

namespace adler {

	namespace detail {
		template <unsigned int base, unsigned int exp>
		struct power {
			enum { value = base*power<base,exp-1>::value };
		};

		template <unsigned int base>
		struct power<base,0> {
			enum { value = 1 };
		};
	
		template <unsigned int dimension, typename Real_>
		struct grad {
			inline Real_ operator () (int hash, Real_ v[dimension])
			{
				BOOST_STATIC_ASSERT(dimension > 3); // specialized for simple dimensions
				unsigned int h = hash % (dimension*power<2,dimension-1>::value);

				unsigned int zeroDimension = h%dimension;

				Real_ ret = 0;

				unsigned int i = 0;

				while (i<zeroDimension)
				{
					ret += ((h&(1<<i)) ==0)?v[i]:-v[i];
					++i;
				}				

				++i; // skip the zero dimension

				while (i<dimension)
				{
					ret += ((h&(1<<i))==0)?v[i]:-v[i];
					++i;
				}

				return ret;

			}
		};

		template <typename Real_>
		struct grad<1,Real_> {
			inline Real_ operator() (int hash, Real_ x)
			{
				int h = hash & 1;
				return (h==0)?x:-x;
			}
		};

		template <typename Real_>
		struct grad<2,Real_> {
			inline Real_ operator() (int hash, Real_ x, Real_ y)
			{
				int h = hash & 3;
				return ((h&1)==0 ? x : -x) + ((h&2)==0?y:-y);
			}
		};

		template <typename Real_>
		struct grad<3,Real_> {
			inline Real_ operator() (int hash, Real_ x, Real_ y, Real_ z)
			{
				int h = hash & 15;
				Real_ u = h<8 ? x : y;
				Real_ v = h<4 ? y : h==12||h==14 ? x : z;
				return ((h&1)==0 ? u : -u) + ((h&2)==0?v:-v);
			}
		};
		
		template <typename Real_>
		inline Real_ fade(Real_ t)
		{
			return t*t*t*(t*(t*6-15)+10);
		}

		template <typename Real_>
		inline Real_ lerp(Real_ t, Real_ a, Real_ b)
		{
			return a + t*(b-a);
		}

		template <unsigned int dimension,
			unsigned int treeheight,
			typename Real_,
			size_t PTableSize_>
		struct lerp_builder {
			inline Real_ operator() (Real_ fades[dimension], Real_ vrel[dimension],
				int floors[dimension], const int p[PTableSize_], unsigned int bitflags)				
			{
				BOOST_STATIC_ASSERT(treeheight < dimension);
				BOOST_STATIC_ASSERT(dimension > 3); // special case for simple dimensions
				
				// if treeheight+1=dimension (i.e. initial call), -> bitflags = 0
				assert( !(treeheight+1==dimension) || (bitflags==0) );
			
				Real_ val1, val2;

				lerp_builder<dimension,treeheight-1,Real_,PTableSize_> l;

				/*for (int i=0; i<dimension-treeheight; ++i) std::clog << '\t';
				std::clog << "lerp(fades[" << treeheight << "],\n";*/

				val1 = l(fades,vrel,floors,p,bitflags);
				vrel[treeheight]--;
				bitflags |= (1<<treeheight);
				
				val2 = l(fades,vrel,floors,p,bitflags);
				vrel[treeheight]++;
				
				/*for (int i=0; i<dimension-treeheight; ++i) std::clog << '\t';
				std::clog << ")\n";*/

				return lerp(fades[treeheight],val1,val2);
			}
		};

		template <unsigned int dimension,
			typename Real_,
			size_t PTableSize_>
		struct lerp_builder<dimension,0,Real_,PTableSize_> {
			inline Real_ operator() (Real_ fades[dimension], Real_ vrel[dimension],
				int floors[dimension], const int p[PTableSize_], unsigned int bitflags)
			{
				// TODO: specialize most inward case, which actually calls the gradients
				grad<dimension, Real_> g;

				Real_ val1, val2;

				int hash = 0;

				/*for (int i=0; i<dimension; ++i) std::clog << '\t';
				std::clog << "lerp(fades[" << 0 << "],\n";*/

				for (unsigned i=0; i<dimension; ++i)
				{
					hash = p[floors[i]+((((1<<i)&bitflags)==0)?0:1)+hash];
				}
				
				val1 = g(hash, vrel);
				vrel[0]--;

				hash = 0;
				bitflags |= 1;
				for (unsigned i=0; i<dimension; ++i)
				{
					hash = p[floors[i]+((((1<<i)&bitflags)==0)?0:1)+hash];
				}
				
				val2 = g(hash, vrel);
				vrel[0]++;				
				
				/*for (int i=0; i<dimension; ++i) std::clog << '\t';
				std::clog << ")\n";*/
				
				return lerp(fades[0],val1,val2);
			}
		};

		
	}

	template <typename RNG_, typename Real_=float, size_t PTableSize_=512>
	class Perlin {
	private:
		int p[PTableSize_];

	public:
		Perlin(RNG_ &rng) {
			BOOST_STATIC_ASSERT( (PTableSize_&1) == 0);

			boost::uniform_int<> bytesize(0,255);
			boost::variate_generator<RNG_,boost::uniform_int<> >
				die(rng,bytesize);
			
			for (unsigned int i=0; i<PTableSize_/2; ++i)
			{
				p[i] = p[i+PTableSize_/2] = die();				
			}
		}

		template <unsigned int dimension>
		Real_ noise(Real_ v[dimension])
		{
			BOOST_STATIC_ASSERT(dimension > 0);

			int floors[dimension];
			Real_ vrel[dimension];

			for (unsigned int i=0; i<dimension; ++i)
			{
				floors[i] = static_cast<int>( std::floor(v[i]) );
				vrel[i] = v[i] - floors[i];
			}

			Real_ fades[dimension];
			for (unsigned int i=0; i<dimension; ++i)
			{
				fades[i] = detail::fade(vrel[i]);
			}


			detail::lerp_builder<dimension,dimension-1,Real_,PTableSize_> l;

			return l(fades, vrel, floors, p, 0);

		}

		template <unsigned int dimension>
		Real_ noise(Real_ x)
		{
			BOOST_STATIC_ASSERT(dimension==1);

			int X = static_cast<int>(std::floor(x)) % (PTableSize_/2);
			x -= std::floor(x);

			Real_ u = detail::fade(x);

			int A=p[X], B=p[X+1];

			detail::grad<1,Real_> g;

			return detail::lerp(u,g(A,x),g(B,x-1));
		}

		template <unsigned int dimension>
		Real_ noise(Real_ x, Real_ y, Real_ z)
		{
			BOOST_STATIC_ASSERT(dimension==3);

			int X = static_cast<int>(std::floor(x)) % (PTableSize_/2);
			int Y = static_cast<int>(std::floor(y)) % (PTableSize_/2);
			int Z = static_cast<int>(std::floor(z)) % (PTableSize_/2);

			x -= std::floor(x);
			y -= std::floor(y);
			z -= std::floor(z);

			Real_ u = detail::fade(x), v = detail::fade(y), w = detail::fade(z);

			int A=p[X]+Y,   AA=p[A]+Z, AB=p[A+1]+Z,
			    B=p[X+1]+Y, BA=p[B]+Z, BB=p[B+1]+Z;

			detail::grad<3,Real_> g;

			return
				detail::lerp(w,
					detail::lerp(v,
						detail::lerp(u, g(p[AA], x,y,z), g(p[BA], x-1,y,z)),
						detail::lerp(u, g(p[AB], x,y-1,z), g(p[BB], x-1,y-1,z))),
					detail::lerp(v,
						detail::lerp(u, g(p[AA+1], x,y,z-1), g(p[BA+1], x-1,y,z-1)),
						detail::lerp(u, g(p[AB+1], x,y-1,z-1), g(p[BB+1], x-1,y-1,z-1))));
		}
		
	};

} // end namespace adler

#endif

