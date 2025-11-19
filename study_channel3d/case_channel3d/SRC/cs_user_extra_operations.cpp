/*============================================================================
 * User extra operations:
 *  - Integrated forces on pillar_walls (rectangular pillar)
 *  - Drag and lift coefficients (Cd, Cl)
 *  - Approximate pressure loss between inlet and outlet (dp)
 *  - Output as CSV time series
 *============================================================================*/

#include "base/cs_defs.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "cs_headers.h"

BEGIN_C_DECLS

/*---------------------------------------------------------------------------*/
/* Time-plot handle                                                           */
/*---------------------------------------------------------------------------*/

static cs_time_plot_t *_pillar_force_plot = NULL;

/*---------------------------------------------------------------------------*/
/* Reference values for coefficients                                         */
/*---------------------------------------------------------------------------*/

static const double rho_ref = 1.0; /* density [kg/m3] */
static const double U_ref = 1.0; /* reference velocity [m/s] (Re = 100) */
static const double D_ref = 1.0; /* characteristic length: pillar height in y [m] */
static const double L_ref = 0.41; /* spanwise thickness in z [m] */
static const double A_ref = D_ref * L_ref; /* reference area [m2] */

/*---------------------------------------------------------------------------*/
/*!
 * \brief Extra operations called at the end of each time step.
 */
/*---------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t  *domain)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  /*--------------------------------------------------------------------*/
  /* 1) Initialize CSV time-plot (once)                                 */
  /*--------------------------------------------------------------------*/

  const int n_variables = 6; /* Fx, Fy, Fz, Cd, Cl, dp */

  if (cs_glob_rank_id < 1 && _pillar_force_plot == NULL) {

    int plot_buffer_steps = -1; /* no limit */
    double plot_flush_wtime = 60.0; /* flush every 60s */
    cs_time_plot_format_t plot_format = CS_TIME_PLOT_CSV;
    bool use_iteration = (ts->is_local) ? true : false;

    const char **labels;
    BFT_MALLOC(labels, n_variables, const char *);

    labels[0] = "Fx";
    labels[1] = "Fy";
    labels[2] = "Fz";
    labels[3] = "Cd";
    labels[4] = "Cl";
    labels[5] = "dp";

    _pillar_force_plot = cs_time_plot_init_probe(
        "pillar_forces",
        "",
        plot_format,
        use_iteration,
        plot_flush_wtime,
        plot_buffer_steps,
        n_variables,
        NULL,
        NULL,
        labels
    );
    BFT_FREE(labels);
  }

  /*--------------------------------------------------------------------*/
  /* 2) Compute total force on pillar_walls                              */
  /*--------------------------------------------------------------------*/

  cs_real_3_t total_force = {0., 0., 0.};

  /* boundary_stress field: traction vector (pressure + viscous) */
  cs_field_t *b_forces = cs_field_by_name_try("boundary_stress");

  if (b_forces != NULL) {

    const cs_real_3_t *bpro_forces =
      (const cs_real_3_t *)(b_forces->val);

    const cs_real_t *b_face_surf =
      domain->mesh_quantities->b_face_surf;

    /* Physical Surface("pillar_walls") in .geo file */
    const cs_zone_t *zw = cs_boundary_zone_by_name("pillar_walls");

    if (zw != NULL) {

      for (cs_lnum_t e_id = 0; e_id < zw->n_elts; e_id++) {

        const cs_lnum_t face_id = zw->elt_ids[e_id];

        /* Force = traction * area */
        for (int i = 0; i < 3; i++)
          total_force[i] += bpro_forces[face_id][i] * b_face_surf[face_id];
      }

      /* Parallel sum */
      cs_parall_sum(3, CS_REAL_TYPE, total_force);
    }
  }

  /*--------------------------------------------------------------------*/
  /* 3) Compute drag and lift coefficients */
  /*    x-direction: drag, y-direction: lift */
  /*--------------------------------------------------------------------*/

  double Cd = 0.0;
  double Cl = 0.0;

  if (rho_ref > 0.0 && U_ref > 0.0 && A_ref > 0.0) {

    const double q_ref = 0.5 * rho_ref * U_ref * U_ref;

    /* Drag: Fx, Lift: Fy */
    Cd = total_force[0] / (q_ref * A_ref);
    Cl = total_force[1] / (q_ref * A_ref);
  }

  /*--------------------------------------------------------------------*/
  /* 4) Approximate pressure loss between inlet and outlet               */
  /*    using normal component of boundary_stress                        */
  /*--------------------------------------------------------------------*/

  double dp = 0.0;

  if (b_forces != NULL) {

    const cs_real_3_t *bpro_forces =
      (const cs_real_3_t *)(b_forces->val);

    const cs_real_t *b_face_surf =
      domain->mesh_quantities->b_face_surf;

    const cs_real_3_t *b_face_normal =
      (const cs_real_3_t *)(domain->mesh_quantities->b_face_normal);

    double sum_tn_area_in  = 0.0, sum_area_in  = 0.0;
    double sum_tn_area_out = 0.0, sum_area_out = 0.0;

    const cs_zone_t *z_in  = cs_boundary_zone_by_name("inlet");
    const cs_zone_t *z_out = cs_boundary_zone_by_name("outlet");

    if (z_in != NULL) {
      for (cs_lnum_t e_id = 0; e_id < z_in->n_elts; e_id++) {

        const cs_lnum_t face_id = z_in->elt_ids[e_id];
        const double area = b_face_surf[face_id];

        if (area <= 0.0)
          continue;

        cs_real_3_t t = {
            bpro_forces[face_id][0],
            bpro_forces[face_id][1],
            bpro_forces[face_id][2]
        };

        cs_real_3_t n = {
            b_face_normal[face_id][0],
            b_face_normal[face_id][1],
            b_face_normal[face_id][2]
        };

        const double nx = n[0] / area;
        const double ny = n[1] / area;
        const double nz = n[2] / area;

        const double tn = t[0]*nx + t[1]*ny + t[2]*nz; /* traction · n̂ */

        sum_tn_area_in += tn * area;
        sum_area_in    += area;
      }
    }

    if (z_out != NULL) {
      for (cs_lnum_t e_id = 0; e_id < z_out->n_elts; e_id++) {

        const cs_lnum_t face_id = z_out->elt_ids[e_id];
        const double area = b_face_surf[face_id];

        if (area <= 0.0)
          continue;

        cs_real_3_t t = {
            bpro_forces[face_id][0],
            bpro_forces[face_id][1],
            bpro_forces[face_id][2]
        };

        cs_real_3_t n = {
            b_face_normal[face_id][0],
            b_face_normal[face_id][1],
            b_face_normal[face_id][2]
        };

        const double nx = n[0] / area;
        const double ny = n[1] / area;
        const double nz = n[2] / area;

        const double tn = t[0]*nx + t[1]*ny + t[2]*nz; /* traction · n̂ */

        sum_tn_area_out += tn * area;
        sum_area_out    += area;
      }
    }

    cs_parall_sum(1, CS_REAL_TYPE, &sum_tn_area_in);
    cs_parall_sum(1, CS_REAL_TYPE, &sum_area_in);
    cs_parall_sum(1, CS_REAL_TYPE, &sum_tn_area_out);
    cs_parall_sum(1, CS_REAL_TYPE, &sum_area_out);

    double t_n_in  = 0.0;
    double t_n_out = 0.0;

    if (sum_area_in  > 0.0)
      t_n_in  = sum_tn_area_in  / sum_area_in;

    if (sum_area_out > 0.0)
      t_n_out = sum_tn_area_out / sum_area_out;

    dp = -(t_n_in - t_n_out);
  }

  /*--------------------------------------------------------------------*/
  /* 5) Write time series to CSV via time-plot                           */
  /*--------------------------------------------------------------------*/

  if (_pillar_force_plot != NULL) {

    double vals[n_variables];

    vals[0] = total_force[0]; /* Fx */
    vals[1] = total_force[1]; /* Fy */
    vals[2] = total_force[2]; /* Fz */
    vals[3] = Cd;
    vals[4] = Cl;
    vals[5] = dp;

    cs_time_plot_vals_write(
        _pillar_force_plot,
        ts->nt_cur,
        ts->t_cur,
        n_variables,
        vals
    );
  }
}

END_C_DECLS
