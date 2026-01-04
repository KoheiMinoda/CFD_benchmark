#include "base/cs_defs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "cs_headers.h"
#include "base/cs_field.h"
#include "mesh/cs_mesh_quantities.h"

BEGIN_C_DECLS

/* -----------------------------------------------------------------------
 * Static plot handles
 * ----------------------------------------------------------------------- */
static cs_time_plot_t *_cylinder_drag_plot = NULL;

/* -----------------------------------------------------------------------
 * Parameters
 * ----------------------------------------------------------------------- */
static const double rho_ref   = 1.0;  // [kg/m^3]
static const double U_ref     = 1.0;  // [m/s]
static const double D_ref     = 1.0;  // [m] Cylinder Diameter
static const double L_span    = 1.0;  // [m] Domain Thickness

/* -----------------------------------------------------------------------
 * Main function
 * ----------------------------------------------------------------------- */
void
cs_user_extra_operations(cs_domain_t *domain)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  /* =================================================================== */
  /* ============= Drag Calculation (Projected Stress Method) ========== */
  /* =================================================================== */
  
  if (cs_glob_rank_id < 1 && _cylinder_drag_plot == NULL) {
    const char *labels[] = {
      "Re", 
      "Cd_pressure", 
      "Cd_viscous",
      "Cd_total", 
      "Cl_total"
    };

    _cylinder_drag_plot = cs_time_plot_init_probe(
      "drag_monitoring", 
      "DragDecomposition", 
      CS_TIME_PLOT_CSV, 
      (ts->is_local), 
      10.0, 
      -1, 
      5, 
      NULL, 
      NULL, 
      labels
    );
  }

  /* Retrieve Boundary Stress Field (Must be enabled in GUI) */
  cs_field_t *f_stress = cs_field_by_name_try("boundary_stress");

  if (f_stress != NULL) {

    /* Pointers to Geometry and Stress */
    /* Note: Casting to handle C++ strict type checking */
    const cs_real_3_t *b_stress_val = (const cs_real_3_t *)(f_stress->val);
    const cs_real_t *b_face_surf    = domain->mesh_quantities->b_face_surf;
    const cs_real_3_t *b_face_normal = (const cs_real_3_t *)domain->mesh_quantities->b_face_normal;

    cs_real_3_t force_total = {0., 0., 0.};
    cs_real_3_t force_pressure = {0., 0., 0.};
    cs_real_3_t force_viscous = {0., 0., 0.};
    
    /* Loop over cylinder walls */
    const cs_zone_t *zn = cs_boundary_zone_by_name("cylinder_walls");
    
    if (zn != NULL) {
      for (cs_lnum_t e_id = 0; e_id < zn->n_elts; e_id++) {
        const cs_lnum_t face_id = zn->elt_ids[e_id];
        
        double area = b_face_surf[face_id];
        
        double nx = b_face_normal[face_id][0] / area;
        double ny = b_face_normal[face_id][1] / area;
        double nz = b_face_normal[face_id][2] / area;

        double fx_total = b_stress_val[face_id][0] * area;
        double fy_total = b_stress_val[face_id][1] * area;
        double fz_total = b_stress_val[face_id][2] * area;

        double f_normal_mag = fx_total*nx + fy_total*ny + fz_total*nz;

        double fx_p = f_normal_mag * nx;
        double fy_p = f_normal_mag * ny;
        double fz_p = f_normal_mag * nz;

        double fx_v = fx_total - fx_p;
        double fy_v = fy_total - fy_p;
        double fz_v = fz_total - fz_p;

        force_total[0] += fx_total;
        force_total[1] += fy_total;
        force_total[2] += fz_total;

        force_pressure[0] += fx_p;
        force_pressure[1] += fy_p;
        force_pressure[2] += fz_p;

        force_viscous[0] += fx_v;
        force_viscous[1] += fy_v;
        force_viscous[2] += fz_v;
      }
    }

    /* Parallel Sum */
    cs_parall_sum(3, CS_REAL_TYPE, force_total);
    cs_parall_sum(3, CS_REAL_TYPE, force_pressure);
    cs_parall_sum(3, CS_REAL_TYPE, force_viscous);

    /* Normalization */
    double q_ref = 0.5 * rho_ref * U_ref * U_ref;
    double area_ref = D_ref * L_span;
    double denom = q_ref * area_ref;
    if (denom < 1e-12) denom = 1.0;

    double Cd_p = force_pressure[0] / denom;
    double Cd_v = force_viscous[0] / denom;
    double Cd_t = force_total[0]    / denom;
    double Cl_t = force_total[1]    / denom;

    /* Get Current Re */
    double current_mu = 0.0;
    cs_field_t *f_mu = cs_field_by_name_try("molecular_viscosity");
    if (f_mu) {
        current_mu = f_mu->val[0]; 
    }
    double current_Re = (current_mu > 1e-9) ? (rho_ref * U_ref * D_ref / current_mu) : 0.0;

    /* Write Output */
    cs_real_t out[5];
    out[0] = current_Re;
    out[1] = Cd_p;
    out[2] = Cd_v;
    out[3] = Cd_t;
    out[4] = Cl_t;

    cs_time_plot_vals_write(
      _cylinder_drag_plot,
      ts->nt_cur, 
      ts->t_cur, 
      5, 
      out
    );
  }
}

END_C_DECLS
