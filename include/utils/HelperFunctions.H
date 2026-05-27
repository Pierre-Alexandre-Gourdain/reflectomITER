#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <string>

#include "utils/Constants.H"

/**
 * @brief Time metadata associated with a solver update.
 *
 * `TimeState` stores the current timestep and the relevant old/current/new
 * times used by time-integration routines or field-coupling logic.
 */
struct TimeState {
    /**
     * @brief Current timestep size.
     */
    amrex::Real dt = 0;

    /**
     * @brief Current simulation time.
     */
    amrex::Real time = 0;

    /**
     * @brief Previous simulation time.
     */
    amrex::Real time_old = 0;

    /**
     * @brief New simulation time after the current update.
     */
    amrex::Real time_new = 0;
};

/**
 * @brief Average an edge-centered field to cell centers.
 *
 * The destination `MultiFab` must already be allocated on a compatible
 * `BoxArray` and `DistributionMapping`.
 *
 * @param mf_edge Source edge-centered field.
 * @param mf_cc Destination cell-centered field.
 * @param edge_dir Edge direction: `0` for x edges, `1` for y edges, `2` for z
 * edges.
 */
void average_edge_to_cellcenter(const amrex::MultiFab& mf_edge,
                                amrex::MultiFab& mf_cc,
                                int edge_dir);

/**
 * @brief Average a cell-centered field to an edge-centered layout.
 *
 * The destination `MultiFab` must already have the nodal staggering appropriate
 * for the requested edge direction.
 *
 * @param mf_cc Source cell-centered field.
 * @param mf_edge Destination edge-centered field.
 * @param edge_dir Edge direction: `0` for x edges, `1` for y edges, `2` for z
 * edges.
 */
void average_cellcenter_to_edge(const amrex::MultiFab& mf_cc,
                                amrex::MultiFab& mf_edge,
                                int edge_dir);

/**
 * @brief Average a face-centered field to cell centers.
 *
 * The destination `MultiFab` must already be allocated on a compatible
 * `BoxArray` and `DistributionMapping`.
 *
 * @param mf_face Source face-centered field.
 * @param mf_cc Destination cell-centered field.
 * @param face_dir Face-normal direction: `0` for x faces, `1` for y faces,
 * `2` for z faces.
 */
void average_face_to_cellcenter(const amrex::MultiFab& mf_face,
                                amrex::MultiFab& mf_cc,
                                int face_dir);

/**
 * @brief Average a cell-centered field to a face-centered layout.
 *
 * The destination `MultiFab` must already have the nodal staggering appropriate
 * for the requested face direction.
 *
 * @param mf_cc Source cell-centered field.
 * @param mf_face Destination face-centered field.
 * @param face_dir Face-normal direction: `0` for x faces, `1` for y faces,
 * `2` for z faces.
 */
void average_cellcenter_to_face(const amrex::MultiFab& mf_cc,
                                amrex::MultiFab& mf_face,
                                int face_dir);

/**
 * @brief Average a node-centered field to cell centers.
 *
 * Each cell-centered value is formed from the eight surrounding nodal values.
 * The destination `MultiFab` must already be allocated.
 *
 * @param mf_node Source node-centered field.
 * @param mf_cc Destination cell-centered field.
 */
void average_node_to_cellcenter(const amrex::MultiFab& mf_node,
                                amrex::MultiFab& mf_cc);

/**
 * @brief Average a cell-centered field to nodes.
 *
 * Each nodal value is formed from the eight surrounding cell-centered values.
 * The destination `MultiFab` must already be allocated.
 *
 * @param mf_cc Source cell-centered field.
 * @param mf_node Destination node-centered field.
 */
void average_cellcenter_to_node(const amrex::MultiFab& mf_cc,
                                amrex::MultiFab& mf_node);

/**
 * @brief Print AMReX box and tile metadata for a field.
 *
 * This diagnostic helper prints valid boxes, tile boxes, grown tile boxes, and
 * FAB allocation boxes. It is useful for debugging ghost-cell access, tiling,
 * and layout mismatches.
 *
 * @param F Field to inspect.
 * @param name Label printed with the field information.
 * @param ngrow Number of ghost cells used when printing grown tile boxes.
 */
void PrintFieldBoxInfo(const amrex::MultiFab& F,
                       const std::string& name,
                       int ngrow = 0);

/**
 * @brief Print overlap information between two `MultiFab` layouts.
 *
 * The routine compares tile boxes from two fields and reports the number of
 * integer grid points in each overlap.
 *
 * @param F1 First field.
 * @param name1 Name printed for the first field.
 * @param F2 Second field.
 * @param name2 Name printed for the second field.
 */
void CheckOverlap(const amrex::MultiFab& F1,
                  const std::string& name1,
                  const amrex::MultiFab& F2,
                  const std::string& name2);

/**
 * @brief Return a geometry with an expanded physical `RealBox`.
 *
 * The index-space domain is preserved, but the physical coordinate bounds are
 * grown by `ngrow` cell widths in each direction.
 *
 * @param geom_in Input geometry.
 * @param ngrow Number of cell widths by which to grow the physical bounds.
 * @return Geometry with enlarged physical bounds.
 */
amrex::Geometry grow_geometry(const amrex::Geometry& geom_in,
                              int ngrow);

/**
 * @brief Split a possibly overlapping `BoxList` into disjoint boxes.
 *
 * The function constructs elementary boxes induced by all input box boundaries
 * and keeps the pieces covered by at least one input box.
 *
 * @param boxes Input box list, possibly overlapping.
 * @return Disjoint box list covering the same union.
 */
amrex::BoxList makeDisjoint(const amrex::BoxList& boxes);

/**
 * @brief Write a `BoxArray` to disk.
 *
 * The output uses AMReX's native stream format and is mainly intended for
 * debugging or reproducing decomposition problems.
 *
 * @param ba BoxArray to write.
 * @param filename Output filename.
 */
void dumpBoxArray(const amrex::BoxArray& ba,
                  const std::string& filename);

/**
 * @brief Write geometry metadata to a human-readable file.
 *
 * The dump includes the index domain, physical `RealBox`, coordinate system,
 * and periodicity flags.
 *
 * @param geom Geometry to write.
 * @param filename Output filename.
 */
void dumpGeometry(const amrex::Geometry& geom,
                  const std::string& filename);

/**
 * @brief Construct a geometry consistent with a new `BoxArray`.
 *
 * The new geometry uses the same cell size, coordinate system, and periodicity
 * as the input geometry. Its physical bounds are shifted to correspond to the
 * new index-space domain.
 *
 * @param geom_in Reference geometry.
 * @param ba_in Original BoxArray associated with `geom_in`.
 * @param ba_new New BoxArray used to define the output domain.
 * @return Geometry compatible with `ba_new`.
 */
amrex::Geometry MakeGeometryFromBoxArray(const amrex::Geometry& geom_in,
                                         const amrex::BoxArray& ba_in,
                                         const amrex::BoxArray& ba_new);

/**
 * @brief Return the largest side length in a `BoxArray`.
 *
 * @param ba BoxArray to inspect.
 * @return Maximum box side length over all boxes and directions.
 */
int getMaxBoxSize(const amrex::BoxArray& ba);

/**
 * @brief Fill non-periodic ghost cells by extrapolation.
 *
 * Periodic directions are skipped. The extrapolation order is selected by
 * `order`:
 *
 * - `0`: constant extrapolation,
 * - `1`: linear extrapolation,
 * - `2` or larger: quadratic extrapolation.
 *
 * @param mf Field whose ghost cells are filled.
 * @param geom Geometry defining the physical domain and periodicity.
 * @param order Extrapolation order.
 */
void FillOpenBoundaries(amrex::MultiFab& mf,
                        const amrex::Geometry& geom,
                        int order = 1);

/**
 * @brief Fill non-periodic ghost cells with a constant Dirichlet value.
 *
 * Only ghost cells outside non-periodic physical boundaries are modified.
 * Periodic directions are skipped.
 *
 * @param mf Field whose ghost cells are filled.
 * @param geom Geometry defining the physical domain and periodicity.
 * @param val Dirichlet value assigned to physical ghost cells.
 */
void FillDirichletBoundaries(amrex::MultiFab& mf,
                             const amrex::Geometry& geom,
                             amrex::Real val);

/**
 * @brief Generate an outer shell around a physical-domain `BoxArray`.
 *
 * The generated shell is intended for PML or other exterior-domain operations.
 * Periodic directions are excluded. The output `BoxArray`, `Geometry`, and
 * `DistributionMapping` are overwritten by the generated shell data.
 *
 * If `pml_thickness > 1`, it is interpreted as a number of cells. If
 * `pml_thickness <= 1`, it is interpreted as a fraction of the domain size.
 *
 * @param ba Input physical-domain BoxArray.
 * @param geom Geometry associated with the physical domain.
 * @param ba_shell Output shell BoxArray.
 * @param geom_shell Output shell geometry.
 * @param dm_shell Output shell distribution mapping.
 * @param pml_thickness Shell thickness in cells or as a domain fraction.
 */
void outer_shell_generation(const amrex::BoxArray& ba,
                            const amrex::Geometry& geom,
                            amrex::BoxArray& ba_shell,
                            amrex::Geometry& geom_shell,
                            amrex::DistributionMapping& dm_shell,
                            amrex::Real pml_thickness = 0.1);

#endif