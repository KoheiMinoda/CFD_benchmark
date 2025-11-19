/* cs_user_initialization.c */

#include "cs_headers.h"

BEGIN_C_DECLS

/*---------------------------------------------------------------------------*/
/* User initialization                                                       */
/*---------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t  *domain)
{
  /* メッシュとメッシュ量 */
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  /* セル中心座標 */
  const cs_real_3_t *cell_cen =
    (const cs_real_3_t *)(mq->cell_cen);

  /* HGN の void_fraction フィールドを取得 */
  cs_field_t *vf = cs_field_by_name_try("void_fraction");
  if (vf == NULL) {
    bft_printf("WARNING: field 'void_fraction' not found.\n");
    return;
  }

  cs_real_t *vf_val = vf->val;

  /* 全セルを走査して y 座標で 0/1 を設定 */
  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    cs_real_t y = cell_cen[c_id][1];  /* [0]=x, [1]=y, [2]=z */

    if (y <= 0.0)
      vf_val[c_id] = 0.0;   /* 水：空気体積分率 0 */
    else
      vf_val[c_id] = 1.0;   /* 空気：空気体積分率 1 */
  }
}

END_C_DECLS
