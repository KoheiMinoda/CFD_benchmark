/* cs_user_source_terms.c : Rayleigh sponge for x > xs */

#include "cs_defs.h"
#include "cs_headers.h"

#include <string.h>

BEGIN_C_DECLS

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  /* 1) velocity 以外の場には何もしない */
  const cs_field_t *f = cs_field_by_id(f_id);
  if (strcmp(f->name, "velocity") != 0)
    return;

  /* 2) メッシュ情報とセル中心座標を取得 */
  const cs_mesh_t            *m  = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_lnum_t    n_cells  = m->n_cells;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)(mq->cell_cen);

  /* 3) st_imp を 3 成分ベクトルとして扱う */
  cs_real_3_t *simp = (cs_real_3_t *)st_imp;

  /* 4) sponge パラメータ (.geo に合わせる) */
  const cs_real_t Lx        = 20.0;  /* タンク全長 (0–Ltot) */
  const cs_real_t xs        = 15.0;  /* sponge 開始位置 x = Lw = 15 */
  const cs_real_t sigma_max = 2.0;   /* 減衰強さ [1/s] → まずは控えめに */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t x = cell_cen[c_id][0];

    if (x > xs) {

      cs_real_t xi = (x - xs) / (Lx - xs);
      if (xi < 0.) xi = 0.;
      if (xi > 1.) xi = 1.;

      /* quadratic ramp */
      const cs_real_t sigma = sigma_max * xi * xi;

      /* du/dt = ... - sigma * u を実現するために
         st_exp = 0, st_imp = -sigma とする */
      simp[c_id][0] -= sigma;
      simp[c_id][1] -= sigma;
      simp[c_id][2] -= sigma;
    }
  }

  /* 今回 st_exp は使わない */
  (void)st_exp;
}

END_C_DECLS
