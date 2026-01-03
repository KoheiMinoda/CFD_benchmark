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
        
        // 1. Get Stress Vector (Force per Unit Area)
        double sigma[3];
        sigma[0] = b_stress_val[face_id][0];
        sigma[1] = b_stress_val[face_id][1];
        sigma[2] = b_stress_val[face_id][2];

        // 2. Get Unit Normal Vector
        double nx = b_face_normal[face_id][0];
        double ny = b_face_normal[face_id][1];
        double nz = b_face_normal[face_id][2];

        // 3. Project Stress onto Normal (Normal Stress = Pressure component)
        // Dot product: S . n
        double sigma_n = sigma[0]*nx + sigma[1]*ny + sigma[2]*nz;

        // Pressure Force Vector (Normal component * Area)
        double f_p[3];
        f_p[0] = sigma_n * nx * area;
        f_p[1] = sigma_n * ny * area;
        f_p[2] = sigma_n * nz * area;

        // Total Force Vector (Stress * Area)
        double f_t[3];
        f_t[0] = sigma[0] * area;
        f_t[1] = sigma[1] * area;
        f_t[2] = sigma[2] * area;

        // Viscous Force Vector (Total - Pressure)
        // This corresponds to the Tangential component
        double f_v[3];
        f_v[0] = f_t[0] - f_p[0];
        f_v[1] = f_t[1] - f_p[1];
        f_v[2] = f_t[2] - f_p[2];

        // Accumulate
        for (int i = 0; i < 3; i++) {
            force_total[i]    += f_t[i];
            force_pressure[i] += f_p[i];
            force_viscous[i]  += f_v[i];
        }
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
