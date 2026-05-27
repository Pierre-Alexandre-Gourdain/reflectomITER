#define _FD_STENCIL(a,ii,jj,kk,dsi) ((a(i+ii,j+jj,k+kk)-a(i,j,k))*(dsi))
#define _VECT_FD_STENCIL(A,ii,jj,kk,dsi,AXIS) ((A(i+ii,j+jj,k+kk,AXIS)-A(i,j,k,AXIS))*(dsi))

#define fwd_dx(a) _FD_STENCIL(a,1,0,0,dxi)
#define fwd_dy(a) _FD_STENCIL(a,0,1,0,dyi)
#define fwd_dz(a) _FD_STENCIL(a,0,0,1,dzi)

#define fwd_DX(A,AXIS) _VECT_FD_STENCIL(A,1,0,0,dxi,AXIS)
#define fwd_DY(A,AXIS) _VECT_FD_STENCIL(A,0,1,0,dyi,AXIS)
#define fwd_DZ(A,AXIS) _VECT_FD_STENCIL(A,0,0,1,dzi,AXIS)

#define bak_dx(a) _FD_STENCIL(a,-1,0,0,-dxi)
#define bak_dy(a) _FD_STENCIL(a,0,-1,0,-dyi)
#define bak_dz(a) _FD_STENCIL(a,0,0,-1,-dzi)

#define bak_DX(A,AXIS) _VECT_FD_STENCIL(A,-1,0,0,-dxi,AXIS)
#define bak_DY(A,AXIS) _VECT_FD_STENCIL(A,0,-1,0,-dyi,AXIS)
#define bak_DZ(A,AXIS) _VECT_FD_STENCIL(A,0,0,-1,-dzi,AXIS)

/* #define d_dx(a,dir) dir##_dx(a)
#define d_dy(a,dir) dir##_dy(a)
#define d_dz(a,dir) dir##_dz(a)

#define D_DX(A,AXIS,dir) dir##_DX(A,AXIS)
#define D_DY(A,AXIS,dir) dir##_DY(A,AXIS)
#define D_DZ(A,AXIS,dir) dir##_DZ(A,AXIS)

#define d2_dx2(a) ((a(i+1,j,k)-2*a(i,j,k)+a(i-1,j,k))*(dxi)*(dxi))
#define d2_dy2(a) ((a(i,j+1,k)-2*a(i,j,k)+a(i,j-1,k))*(dyi)*(dyi))
#define d2_dz2(a) ((a(i,j,k+1)-2*a(i,j,k)+a(i,j,k-1))*(dzi)*(dzi))

#define D2_DX2(A) ((A(i+1,j,k,X)-2*A(i,j,k,X)+A(i-1,j,k,X))*(dxi)*(dxi))
#define D2_DY2(A) ((A(i,j+1,k,Y)-2*A(i,j,k,Y)+A(i,j-1,k,Y))*(dyi)*(dyi))
#define D2_DZ2(A) ((A(i,j,k+1,Z)-2*A(i,j,k,Z)+A(i,j,k-1,Z))*(dzi)*(dzi)) */