#include <utils/HelperFunctions.H>

#include <fstream>

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

using namespace amrex;

void average_edge_to_cellcenter(const MultiFab& mf_edge,
                                MultiFab& mf_cc,
                                int edge_dir)
{
    /*
     * Average an edge-centered field onto cell centers. The destination
     * MultiFab is assumed to be allocated already on a compatible BoxArray and
     * DistributionMapping.
     *
     * edge_dir selects the direction of the edge-centered quantity:
     *   0 -> x-directed edge,
     *   1 -> y-directed edge,
     *   2 -> z-directed edge.
     */
    for (MFIter mfi(mf_edge); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto arr_edge = mf_edge.const_array(mfi);
        auto arr_cc   = mf_cc.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            if (edge_dir == 0) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_edge(i,j,k)
                    + arr_edge(i,j+1,k)
                    + arr_edge(i,j,k+1)
                    + arr_edge(i,j+1,k+1)
                );
            } else if (edge_dir == 1) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_edge(i,j,k)
                    + arr_edge(i+1,j,k)
                    + arr_edge(i,j,k+1)
                    + arr_edge(i+1,j,k+1)
                );
            } else if (edge_dir == 2) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_edge(i,j,k)
                    + arr_edge(i+1,j,k)
                    + arr_edge(i,j+1,k)
                    + arr_edge(i+1,j+1,k)
                );
            }
        });
    }
}

void average_cellcenter_to_edge(const MultiFab& mf_cc,
                                MultiFab& mf_edge,
                                int edge_dir)
{
    /*
     * Average a cell-centered scalar field onto an edge-centered layout. The
     * target MultiFab must already have the nodal staggering appropriate for
     * the requested edge direction.
     */
    for (MFIter mfi(mf_cc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        auto arr_cc   = mf_cc.const_array(mfi);
        auto arr_edge = mf_edge.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            if (edge_dir == 0) {
                arr_edge(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i-1,j,k));
            } else if (edge_dir == 1) {
                arr_edge(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i,j-1,k));
            } else if (edge_dir == 2) {
                arr_edge(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i,j,k-1));
            }
        });
    }
}

void average_face_to_cellcenter(const MultiFab& mf_face,
                                MultiFab& mf_cc,
                                int face_dir)
{
    /*
     * Average a face-centered quantity onto cell centers. The averaging plane
     * is perpendicular to the face-normal direction.
     */
    for (MFIter mfi(mf_face); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto arr_face = mf_face.const_array(mfi);
        auto arr_cc   = mf_cc.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            if (face_dir == 0) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_face(i,j,k)
                    + arr_face(i,j+1,k)
                    + arr_face(i,j,k+1)
                    + arr_face(i,j+1,k+1)
                );
            } else if (face_dir == 1) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_face(i,j,k)
                    + arr_face(i+1,j,k)
                    + arr_face(i,j,k+1)
                    + arr_face(i+1,j,k+1)
                );
            } else if (face_dir == 2) {
                arr_cc(i,j,k) = 0.25 * (
                      arr_face(i,j,k)
                    + arr_face(i+1,j,k)
                    + arr_face(i,j+1,k)
                    + arr_face(i+1,j+1,k)
                );
            }
        });
    }
}

void average_cellcenter_to_face(const MultiFab& mf_cc,
                                MultiFab& mf_face,
                                int face_dir)
{
    /*
     * Average a cell-centered scalar field onto a face-centered layout. The
     * target MultiFab must already use the nodal staggering appropriate for the
     * requested face direction.
     */
    for (MFIter mfi(mf_cc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        auto arr_cc   = mf_cc.const_array(mfi);
        auto arr_face = mf_face.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            if (face_dir == 0) {
                arr_face(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i-1,j,k));
            } else if (face_dir == 1) {
                arr_face(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i,j-1,k));
            } else if (face_dir == 2) {
                arr_face(i,j,k) = 0.5 * (arr_cc(i,j,k) + arr_cc(i,j,k-1));
            }
        });
    }
}

void average_node_to_cellcenter(const MultiFab& mf_node,
                                MultiFab& mf_cc)
{
    /*
     * Average a node-centered MultiFab onto cell centers. Each cell-centered
     * value is the arithmetic average of the eight surrounding node values.
     */
    for (MFIter mfi(mf_node); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto arr_node = mf_node.const_array(mfi);
        auto arr_cc   = mf_cc.array(mfi);
        auto n_comp   = mf_cc.nComp();

        ParallelFor(bx, n_comp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            arr_cc(i,j,k,n) = 0.125 * (
                  arr_node(i  , j  , k  , n)
                + arr_node(i+1, j  , k  , n)
                + arr_node(i  , j+1, k  , n)
                + arr_node(i+1, j+1, k  , n)
                + arr_node(i  , j  , k+1, n)
                + arr_node(i+1, j  , k+1, n)
                + arr_node(i  , j+1, k+1, n)
                + arr_node(i+1, j+1, k+1, n)
            );
        });
    }
}

void average_cellcenter_to_node(const MultiFab& mf_cc,
                                MultiFab& mf_node)
{
    /*
     * Average a cell-centered MultiFab onto nodes. Each node receives the
     * arithmetic average of the eight surrounding cell-centered values.
     */
    for (MFIter mfi(mf_cc); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto arr_cc   = mf_cc.const_array(mfi);
        auto arr_node = mf_node.array(mfi);
        auto n_comp   = mf_cc.nComp();

        ParallelFor(bx, n_comp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            arr_node(i,j,k,n) = 0.125 * (
                  arr_cc(i-1,j-1,k-1,n)
                + arr_cc(i  ,j-1,k-1,n)
                + arr_cc(i-1,j  ,k-1,n)
                + arr_cc(i  ,j  ,k-1,n)
                + arr_cc(i-1,j-1,k  ,n)
                + arr_cc(i  ,j-1,k  ,n)
                + arr_cc(i-1,j  ,k  ,n)
                + arr_cc(i  ,j  ,k  ,n)
            );
        });
    }
}

void PrintFieldBoxInfo(const MultiFab& F,
                       const std::string& name,
                       int ngrow)
{
    /*
     * Diagnostic routine for inspecting AMReX box, tile, grown-tile, and FAB
     * extents. This is useful when debugging ghost-cell access, tiling behavior,
     * or mismatches between valid and allocated regions.
     */
    Print() << "\n========== Field: " << name << " ==========\n";

    for (MFIter mfi(F, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx    = mfi.validbox();
        const Box& tbx    = mfi.tilebox();
        const Box& ebx    = enclosedCells(mfi.tilebox());
        const Box& gbx    = mfi.growntilebox(ngrow);
        const Box& fabBox = F[mfi].box();

        Print() << "Tile " << mfi.tileIndex() << ":\n"
                << "  number of points     : " << vbx.numPts() << "\n"
                << "  validbox             : " << vbx << "\n"
                << "  tilebox              : " << tbx << "\n"
                << "  enclosed             : " << ebx << "\n"
                << "  growntilebox         : " << gbx << "\n"
                << "  fab.box()            : " << fabBox << "\n"
                << "  smallEnd(tilebox)    : " << tbx.smallEnd()
                << "  bigEnd(tilebox)      : " << tbx.bigEnd() << "\n";
    }

    Print() << "===========================================\n";
}

int NumOverlapPoints(const Box& a, const Box& b)
{
    /*
     * Return the number of integer grid points in the intersection of two
     * AMReX boxes. A non-overlapping pair contributes zero points.
     */
    Box intersect_box = a & b;

    if (intersect_box.ok()) {
        return intersect_box.numPts();
    }

    return 0;
}

void CheckOverlap(const MultiFab& F1,
                  const std::string& name1,
                  const MultiFab& F2,
                  const std::string& name2)
{
    /*
     * Report overlaps between the tile boxes of two MultiFabs. The boxes are
     * converted to enclosed cell boxes before comparison so the output is easier
     * to interpret for mixed index types.
     */
    Vector<Box> boxes_F2;

    for (MFIter mfi(F2, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        boxes_F2.push_back(enclosedCells(mfi.tilebox()));
    }

    for (MFIter mfi(F1, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& b1 = enclosedCells(mfi.tilebox());

        for (int j = 0; j < boxes_F2.size(); ++j) {
            const Box& b2 = boxes_F2[j];

            if (b1.intersects(b2)) {
                int npts = NumOverlapPoints(b1, b2);

                Print() << "Overlap between "
                        << name1 << "[" << mfi.index() << "] and "
                        << name2 << "[" << j << "]: "
                        << npts << " points\n";
            }
        }
    }
}

Geometry grow_geometry(const Geometry& geom_in,
                       int ngrow)
{
    /*
     * Grow the physical RealBox associated with a Geometry while preserving its
     * index-space domain. This is useful when the physical coordinate extent
     * must include ghost-like padding without changing the AMReX domain box.
     */
    Geometry geom_out = geom_in;

    RealBox rb = geom_out.ProbDomain();
    const auto dx = geom_out.CellSizeArray();

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        rb.setLo(d, rb.lo(d) - ngrow * dx[d]);
        rb.setHi(d, rb.hi(d) + ngrow * dx[d]);
    }

    geom_out.ProbDomain(rb);

    return geom_out;
}

BoxList makeDisjoint(const BoxList& boxes)
{
    /*
     * Split an overlapping BoxList into a set of non-overlapping boxes by
     * constructing all elementary cuboids induced by the input box boundaries.
     * This is conservative and robust for shell construction, but it can create
     * many boxes if the input list has many distinct boundary coordinates.
     */
    std::vector<int> xs;
    std::vector<int> ys;
    std::vector<int> zs;

    for (const auto& b : boxes) {
        xs.push_back(b.smallEnd(0));
        xs.push_back(b.bigEnd(0) + 1);

        ys.push_back(b.smallEnd(1));
        ys.push_back(b.bigEnd(1) + 1);

        zs.push_back(b.smallEnd(2));
        zs.push_back(b.bigEnd(2) + 1);
    }

    auto uniquify = [](std::vector<int>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    };

    uniquify(xs);
    uniquify(ys);
    uniquify(zs);

    BoxList disjoint;

    for (size_t ix = 0; ix + 1 < xs.size(); ++ix) {
        for (size_t iy = 0; iy + 1 < ys.size(); ++iy) {
            for (size_t iz = 0; iz + 1 < zs.size(); ++iz) {
                int x0 = xs[ix];
                int x1 = xs[ix+1] - 1;

                int y0 = ys[iy];
                int y1 = ys[iy+1] - 1;

                int z0 = zs[iz];
                int z1 = zs[iz+1] - 1;

                /*
                 * Add the elementary cuboid only if it intersects at least one
                 * of the input boxes.
                 */
                for (const auto& b : boxes) {
                    if (x1 < b.smallEnd(0) || x0 > b.bigEnd(0)) continue;
                    if (y1 < b.smallEnd(1) || y0 > b.bigEnd(1)) continue;
                    if (z1 < b.smallEnd(2) || z0 > b.bigEnd(2)) continue;

                    disjoint.push_back(Box(IntVect(x0,y0,z0),
                                           IntVect(x1,y1,z1)));
                    break;
                }
            }
        }
    }

    return disjoint;
}

void dumpBoxArray(const BoxArray& ba,
                  const std::string& filename)
{
    /*
     * Write an AMReX BoxArray to disk in AMReX's native stream format. This is
     * mainly a debugging helper for reproducing domain-decomposition issues.
     */
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        Abort("Cannot open file for BoxArray dump: " + filename);
    }

    ba.writeOn(ofs);
}

void dumpGeometry(const Geometry& geom,
                  const std::string& filename)
{
    /*
     * Write the main Geometry metadata to a human-readable text file. This is
     * useful when checking whether derived geometries preserve the intended
     * domain, RealBox, coordinate system, and periodicity.
     */
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        Abort("Cannot open file for Geometry dump: " + filename);
    }

    const Box& domain = geom.Domain();
    const RealBox& rb = geom.ProbDomain();
    const int coord_sys = geom.Coord();

    ofs << "Geometry dump\n";
    ofs << "Domain: " << domain << "\n";
    ofs << "RealBox: [";

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        ofs << rb.lo(d) << ", " << rb.hi(d);
        if (d < AMREX_SPACEDIM - 1) ofs << "; ";
    }

    ofs << "]\n";

    ofs << "CoordSys: " << coord_sys << "\n";
    ofs << "Periodic: ";

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        ofs << geom.isPeriodic(d) << " ";
    }

    ofs << "\n";
}

Geometry MakeGeometryFromBoxArray(const Geometry& geom_in,
                                  const BoxArray& ba_in,
                                  const BoxArray& ba_new)
{
    /*
     * Construct a Geometry for a new BoxArray while preserving the cell size,
     * coordinate system, and periodicity of an existing Geometry. The physical
     * RealBox is shifted so the new index-space domain maps consistently into
     * the original coordinate system.
     */
    Box domain_new = enclosedCells(ba_new.minimalBox());
    const Box& domain_in = geom_in.Domain();

    const auto dx = geom_in.CellSizeArray();
    const auto prob_lo_in = geom_in.ProbLoArray();

    RealBox rb_new;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Real lo = prob_lo_in[idim]
                + (domain_new.smallEnd(idim) - domain_in.smallEnd(idim))
                * dx[idim];

        Real hi = lo + domain_new.length(idim) * dx[idim];

        rb_new.setLo(idim, lo);
        rb_new.setHi(idim, hi);
    }

    const int coord = geom_in.Coord();
    const auto is_per = geom_in.isPeriodic();

    Geometry geom_new(domain_new, &rb_new, coord, is_per.data());

    return geom_new;
}

int getMaxBoxSize(const BoxArray& ba)
{
    /*
     * Return the largest side length among all boxes in a BoxArray. This is
     * used to preserve a comparable box-size limit when constructing derived
     * BoxArrays such as PML shells.
     */
    int max_side = 0;

    for (int i = 0; i < ba.size(); ++i) {
        const Box& b = ba[i];
        const IntVect len = b.length();

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            max_side = std::max(max_side, len[d]);
        }
    }

    return max_side;
}

void outer_shell_generation(const BoxArray& ba,
                            const Geometry& geom,
                            BoxArray& ba_shell,
                            Geometry& geom_shell,
                            DistributionMapping& dm_shell,
                            amrex::Real pml_thickness)
{
    /*
     * Build an outer shell BoxArray around a physical-domain BoxArray. The shell
     * is used for PML or other exterior-domain operations. Periodic directions
     * are excluded because no absorbing layer should be generated there.
     */
    using namespace amrex;

    bool is_periodic_x = geom.isPeriodic(0);
    bool is_periodic_y = geom.isPeriodic(1);
    bool is_periodic_z = geom.isPeriodic(2);

    const Box domain = ba.minimalBox();

    IntVect thickness(0);

    if (pml_thickness > 1.0) {
        /*
         * Values larger than one are interpreted as an explicit number of PML
         * cells.
         */
        const int n_pml = static_cast<int>(std::round(pml_thickness));

        thickness[0] = n_pml;
        thickness[1] = n_pml;
        thickness[2] = n_pml;
    } else {
        /*
         * Fractional values are interpreted relative to the domain length. A
         * minimum fraction prevents an accidentally tiny PML on large domains.
         */
        const Real frac = std::max(std::abs(pml_thickness), Real(0.05));

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (!geom.isPeriodic(d) && domain.length(d) > 1) {
                thickness[d] = static_cast<int>(
                    std::ceil(frac * static_cast<Real>(domain.length(d)))
                );

                /*
                 * Current implementation uses a fixed 30-cell shell in each
                 * non-periodic active direction. Remove this assignment if the
                 * fractional thickness above should be used directly.
                 */
                thickness[d] = 30;
            }

            Print() << thickness[d] << std::endl;
        }
    }

    if (is_periodic_x) thickness[0] = 0;
    if (is_periodic_y) thickness[1] = 0;
    if (is_periodic_z) thickness[2] = 0;

    if (is_periodic_x && is_periodic_y && is_periodic_z) {
        thickness = IntVect(0);
    }

    /*
     * Default outputs are the original domain. They are overwritten below only
     * if a nonzero shell is required.
     */
    ba_shell = ba;
    geom_shell = geom;
    dm_shell = DistributionMapping(ba_shell);

    if (thickness[0] == 0 && thickness[1] == 0 && thickness[2] == 0) {
        return;
    }

    BoxList bl_shell;

    Box extended_x = domain;
    if (!is_periodic_x) extended_x.grow(0, thickness[0]);

    Box extended_y = domain;
    if (!is_periodic_y) extended_y.grow(1, thickness[1]);

    Box extended_z = domain;
    if (!is_periodic_z) extended_z.grow(2, thickness[2]);

    /*
     * Create slabs adjacent to each non-periodic side. Each slab is grown in
     * the transverse non-periodic directions so corners and edges are included.
     */
    if (!is_periodic_x) {
        Box extended_yz = extended_y;
        if (!is_periodic_z) extended_yz.grow(2, thickness[2]);

        bl_shell.push_back(adjCellHi(extended_yz, 0, thickness[0]));
        bl_shell.push_back(adjCellLo(extended_yz, 0, thickness[0]));
    }

    if (!is_periodic_y) {
        Box extended_zx = extended_z;
        if (!is_periodic_x) extended_zx.grow(0, thickness[0]);

        bl_shell.push_back(adjCellHi(extended_zx, 1, thickness[1]));
        bl_shell.push_back(adjCellLo(extended_zx, 1, thickness[1]));
    }

    if (!is_periodic_z) {
        Box extended_xy = extended_x;
        if (!is_periodic_y) extended_xy.grow(1, thickness[1]);

        bl_shell.push_back(adjCellHi(extended_xy, 2, thickness[2]));
        bl_shell.push_back(adjCellLo(extended_xy, 2, thickness[2]));
    }

    /*
     * The slabs overlap at edges and corners. Split them into a disjoint list so
     * the resulting BoxArray has unambiguous ownership.
     */
    bl_shell = makeDisjoint(bl_shell);

    ba_shell.define(bl_shell);
    ba_shell.maxSize(getMaxBoxSize(ba));

    dm_shell = DistributionMapping(ba_shell);
    geom_shell = MakeGeometryFromBoxArray(geom, ba, ba_shell);
}

void FillOpenBoundaries(MultiFab& mf,
                        const Geometry& geom,
                        int order)
{
    /*
     * Fill non-periodic ghost cells by extrapolating from the nearest interior
     * cells. The extrapolation order is:
     *
     *   order = 0 -> constant extrapolation,
     *   order = 1 -> linear extrapolation,
     *   order >= 2 -> quadratic extrapolation.
     *
     * Periodic directions are skipped because AMReX periodic FillBoundary
     * should handle those.
     */
    const Box& domain = geom.Domain();
    const IntVect ng  = mf.nGrowVect();

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        auto arr = fab.array();

        const Box& fabbox = fab.box();
        const Box& vbx    = mfi.validbox();

        const auto vlo = lbound(vbx);
        const auto vhi = ubound(vbx);

        const int nc = mf.nComp();

        auto fill_face = [&] (int dir, bool is_lo)
        {
            if (geom.isPeriodic(dir)) return;

            int dom = is_lo ? domain.smallEnd(dir) : domain.bigEnd(dir);
            int sign = is_lo ? +1 : -1;

            IntVect lo(vlo.x, vlo.y, vlo.z);
            IntVect hi(vhi.x, vhi.y, vhi.z);

            lo[dir] = is_lo ? dom - ng[dir] : dom + 1;
            hi[dir] = is_lo ? dom - 1       : dom + ng[dir];

            Box gbx(lo, hi);
            gbx &= fabbox;

            if (!gbx.ok()) return;

            ParallelFor(gbx, nc,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                int idx[3] = {i, j, k};
                int s = is_lo ? (dom - idx[dir]) : (idx[dir] - dom);

                idx[dir] = dom;
                Real u0 = arr(idx[0], idx[1], idx[2], n);

                if (order == 0) {
                    arr(i,j,k,n) = u0;
                    return;
                }

                idx[dir] = dom + sign;
                Real u1 = arr(idx[0], idx[1], idx[2], n);

                if (order == 1) {
                    arr(i,j,k,n) = u0 + s * (u0 - u1);
                    return;
                }

                idx[dir] = dom + 2 * sign;
                Real u2 = arr(idx[0], idx[1], idx[2], n);

                /*
                 * Quadratic extrapolation from the boundary cell and the first
                 * two interior cells.
                 */
                Real a =  0.5_rt * (s + 1) * (s + 2);
                Real b = -1.0_rt * s * (s + 2);
                Real c =  0.5_rt * s * (s + 1);

                arr(i,j,k,n) = a * u0 + b * u1 + c * u2;
            });
        };

        fill_face(0, true);
        fill_face(0, false);

#if AMREX_SPACEDIM >= 2
        fill_face(1, true);
        fill_face(1, false);
#endif

#if AMREX_SPACEDIM == 3
        fill_face(2, true);
        fill_face(2, false);
#endif
    }
}

void FillDirichletBoundaries(MultiFab& mf,
                             const Geometry& geom,
                             Real val)
{
    /*
     * Fill ghost cells outside non-periodic domain faces with a constant
     * Dirichlet value. Only physical boundary ghost cells are touched; periodic
     * directions are skipped.
     */
    const Box& domain = geom.Domain();
    const IntVect ng  = mf.nGrowVect();

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        auto arr = fab.array();

        const Box& fabbox = fab.box();
        const Box& vbx    = mfi.validbox();

        int nc = mf.nComp();

        const auto vlo = lbound(vbx);
        const auto vhi = ubound(vbx);

        if (!geom.isPeriodic(0)) {
            int dom_lo = domain.smallEnd(0);

            Box xlo(Box(IntVect(dom_lo - ng[0], vlo.y, vlo.z),
                        IntVect(dom_lo - 1,     vhi.y, vhi.z)));

            xlo &= fabbox;

            if (xlo.ok()) {
                ParallelFor(xlo, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }

            int dom_hi = domain.bigEnd(0);

            Box xhi(Box(IntVect(dom_hi + 1,     vlo.y, vlo.z),
                        IntVect(dom_hi + ng[0], vhi.y, vhi.z)));

            xhi &= fabbox;

            if (xhi.ok()) {
                ParallelFor(xhi, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }
        }

#if AMREX_SPACEDIM >= 2
        if (!geom.isPeriodic(1)) {
            int dom_lo = domain.smallEnd(1);
            int dom_hi = domain.bigEnd(1);

            Box ylo(Box(IntVect(vlo.x, dom_lo - ng[1], vlo.z),
                        IntVect(vhi.x, dom_lo - 1,     vhi.z)));

            ylo &= fabbox;

            if (ylo.ok()) {
                ParallelFor(ylo, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }

            Box yhi(Box(IntVect(vlo.x, dom_hi + 1,     vlo.z),
                        IntVect(vhi.x, dom_hi + ng[1], vhi.z)));

            yhi &= fabbox;

            if (yhi.ok()) {
                ParallelFor(yhi, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }
        }
#endif

#if AMREX_SPACEDIM == 3
        if (!geom.isPeriodic(2)) {
            int dom_lo = domain.smallEnd(2);
            int dom_hi = domain.bigEnd(2);

            Box zlo(Box(IntVect(vlo.x, vlo.y, dom_lo - ng[2]),
                        IntVect(vhi.x, vhi.y, dom_lo - 1)));

            zlo &= fabbox;

            if (zlo.ok()) {
                ParallelFor(zlo, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }

            Box zhi(Box(IntVect(vlo.x, vlo.y, dom_hi + 1),
                        IntVect(vhi.x, vhi.y, dom_hi + ng[2])));

            zhi &= fabbox;

            if (zhi.ok()) {
                ParallelFor(zhi, nc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    arr(i,j,k,n) = val;
                });
            }
        }
#endif
    }
}