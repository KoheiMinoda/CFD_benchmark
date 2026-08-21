#include "base/cs_defs.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "cs_headers.h"

BEGIN_C_DECLS

static cs_time_plot_t *_force_plot = NULL;

/* ---- Case 043 reference values (real scale, water) ---- */
static const double rho_ref = 998.0;      /* [kg/m^3] */
static const double U_ref   = 0.535;      /* [m/s] inlet velocity */
static const double H_ref   = 0.04;       /* [m] cylinder side */
static const double Lz_ref  = M_PI*0.04;  /* [m] spanwise length (pi*H) */

/* projected (frontal) area for CD/CL normalization:
   height H  x  spanwise length Lz                                  */
static const double A_ref = 0.04 * (M_PI*0.04);   /* H * Lz */

void
cs_user_extra_operations(cs_domain_t *domain)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  /* ---- init CSV writer (rank 0) ---- */
  if (cs_glob_rank_id == 0 && _force_plot == NULL) {

    const char *labels[] = {
      "Fx", "Fy", "Fz",          /* dimensional forces [N] */
      "Cd", "Cl", "Cz",          /* force coefficients    */
      "yplus_max", "yplus_mean"  /* wall resolution monitor */
    };

    _force_plot = cs_time_plot_init_probe(
      "force_tracking", "force_tracking",
      CS_TIME_PLOT_CSV,
      0, 10.0, -1,
      8,                 /* number of columns */
      NULL, NULL, labels);
  }

  /* ================================================================
     Force integration on cylinder walls
     ================================================================ */
  {
    cs_field_t *b_forces = cs_field_by_name_try("boundary_stress");
    cs_field_t *f_yplus  = cs_field_by_name_try("yplus");

    cs_real_3_t total_force = {0., 0., 0.};
    double yplus_max_loc = 0.0;
    double yplus_sum_loc = 0.0;
    double area_sum_loc  = 0.0;

    if (b_forces != NULL) {

      const cs_real_t   *b_face_surf =
        domain->mesh_quantities->b_face_surf;
      const cs_real_3_t *bpro_forces =
        (const cs_real_3_t *)(b_forces->val);

      const char *target_zones[] = {
        "cyl_top",
        "cyl_bottom",
        "cyl_limits_xmin",
        "cyl_limits_xmax"
      };
      const int num_zones = 4;

      for (int z = 0; z < num_zones; z++) {
        const cs_zone_t *zn = cs_boundary_zone_by_name(target_zones[z]);
        if (zn == NULL) continue;

        for (cs_lnum_t e_id = 0; e_id < zn->n_elts; e_id++) {
          const cs_lnum_t face_id = zn->elt_ids[e_id];
          const double Sf = b_face_surf[face_id];

          /* boundary_stress is force per unit area -> multiply by area */
          for (int i = 0; i < 3; i++)
            total_force[i] += bpro_forces[face_id][i] * Sf;

          const double yplus =
            (f_yplus != NULL) ? f_yplus->val[face_id] : 0.0;
          if (yplus > yplus_max_loc) yplus_max_loc = yplus;
          yplus_sum_loc += yplus * Sf;
          area_sum_loc  += Sf;
        }
      }
    }

    /* parallel reduction */
    cs_parall_sum(3, CS_REAL_TYPE, total_force);
    cs_parall_max(1, CS_DOUBLE, &yplus_max_loc);

    double yp_budget[2] = {yplus_sum_loc, area_sum_loc};
    cs_parall_sum(2, CS_DOUBLE, yp_budget);

    const double yplus_max = yplus_max_loc;
    const double yplus_mean =
      (yp_budget[1] > 1e-30) ? yp_budget[0]/yp_budget[1] : 0.0;

    /* ---- write to CSV (rank 0) ---- */
    if (cs_glob_rank_id == 0 && b_forces != NULL) {

      const double q = 0.5 * rho_ref * U_ref * U_ref * A_ref; /* 1/2 rho U^2 A */

      cs_real_t out[8] = {0.};
      out[0] = total_force[0];
      out[1] = total_force[1];
      out[2] = total_force[2];
      if (q > 1e-30) {
        out[3] = total_force[0] / q;   /* Cd */
        out[4] = total_force[1] / q;   /* Cl */
        out[5] = total_force[2] / q;   /* Cz */
      }
      out[6] = (cs_real_t)yplus_max;
      out[7] = (cs_real_t)yplus_mean;

      cs_time_plot_vals_write(
        _force_plot, ts->nt_cur, ts->t_cur, 8, out);
    }
  }

  /* ================================================================
     Vorticity field (for POD / wake visualization)
     ================================================================ */
  {
    cs_field_t *vort = cs_field_by_name_try("vorticity");

    if (vort != NULL && vort->dim == 3) {
      const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
      const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

      cs_real_33_t *gradv;
      BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);
      cs_field_gradient_vector(CS_F_(vel), false, 1, gradv);

      cs_real_3_t *omega = (cs_real_3_t *)(vort->val);
      for (cs_lnum_t c = 0; c < n_cells; c++) {
        omega[c][0] = gradv[c][1][2] - gradv[c][2][1];
        omega[c][1] = gradv[c][2][0] - gradv[c][0][2];
        omega[c][2] = gradv[c][0][1] - gradv[c][1][0];
      }

      BFT_FREE(gradv);
    }
  }
}

END_C_DECLS
