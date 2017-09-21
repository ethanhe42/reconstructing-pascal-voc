/* $Id: mesh_run.c,v 1.28 2004/04/30 07:50:21 aspert Exp $ */


/*
 *
 *  Copyright (C) 2001-2004 EPFL (Swiss Federal Institute of Technology,
 *  Lausanne) This program is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA.
 *
 *  In addition, as a special exception, EPFL gives permission to link
 *  the code of this program with the Qt non-commercial edition library
 *  (or with modified versions of Qt non-commercial edition that use the
 *  same license as Qt non-commercial edition), and distribute linked
 *  combinations including the two.  You must obey the GNU General
 *  Public License in all respects for all of the code used other than
 *  Qt non-commercial edition.  If you modify this file, you may extend
 *  this exception to your version of the file, but you are not
 *  obligated to do so.  If you do not wish to do so, delete this
 *  exception statement from your version.
 *
 *  Authors : Nicolas Aspert, Diego Santa-Cruz and Davy Jacquet
 *
 *  Web site : http://mesh.epfl.ch
 *
 *  Reference :
 *   "MESH : Measuring Errors between Surfaces using the Hausdorff distance"
 *   in Proceedings of IEEE Intl. Conf. on Multimedia and Expo (ICME) 2002, 
 *   vol. I, pp. 705-708, available on http://mesh.epfl.ch
 *
 */







#include <time.h>
#include <string.h>
#include <xalloc.h>
#include <model_analysis.h>
#include <compute_error.h>
#include <model_in.h>
#include <geomutils.h>

#include <mesh_run.h>

/* Reads a model from file 'fname' and returns the model read. If an error
 * occurs a message is printed and the program exists. 
 */
static struct model *read_model_file(const char *fname)
{
  int rcode;
  struct model *m;
  const char *errstr;
  
  rcode = read_fmodel(&m,fname,MESH_FF_AUTO,1);
  if (rcode <= 0) {
    switch (rcode) {
    case 0:
      errstr = "no models in file";
      break;
    case MESH_NO_MEM:
      errstr = "no memory";
      break;
    case MESH_CORRUPTED:
      errstr = "corrupted file or I/O error";
      break;
    case MESH_MODEL_ERR:
      errstr = "model error";
      break;
    case MESH_NOT_TRIAG:
      errstr = "not a triangular mesh model";
      break;
    case MESH_BAD_FF:
      errstr = "unrecognized file format";
      break;
    case MESH_BAD_FNAME:
      errstr = strerror(errno);
      break;
    default:
      errstr = "unknown error";
    }
    fprintf(stderr,"ERROR: %s: %s\n",fname,errstr);
    exit(1);
  } else if (m->num_faces == 0) {
    fprintf(stderr,"ERROR: %s: empty model (no faces)\n",fname);
    exit(1);
  }
  return m;
}

/* see mesh_run.h */
double mesh_run(const struct args *args, struct model_error *model1,
              struct model_error *model2, struct outbuf *out,
              struct prog_reporter *progress)
{
  clock_t start_time;
  struct dist_surf_surf_stats stats;
  struct dist_surf_surf_stats stats_rev;
  double bbox1_diag,bbox2_diag;
  struct model_info *m1info,*m2info;
  double abs_sampling_step,abs_sampling_dens;
  int nv_empty,nf_empty;
  double output;
  
  /* Read models from input files */
  memset(model1,0,sizeof(*model1));
  memset(model2,0,sizeof(*model2));
  m1info = (struct model_info*) xa_malloc(sizeof(*m1info));
  m2info = (struct model_info*) xa_malloc(sizeof(*m2info));
  start_time = clock();
  model1->mesh = read_model_file(args->m1_fname);
  start_time = clock();
  model2->mesh = read_model_file(args->m2_fname);

  /* Analyze models (we don't need normals for model 1, so we don't request
   * for it to be oriented). */
  start_time = clock();
  bbox1_diag = dist_v(&model1->mesh->bBox[0], &model1->mesh->bBox[1]);
  bbox2_diag = dist_v(&model2->mesh->bBox[0], &model2->mesh->bBox[1]);
  analyze_model(model1->mesh,m1info,0,args->verb_analysis,out,"model 1");
  model1->info = m1info;
  analyze_model(model2->mesh,m2info,1,args->verb_analysis,out,"model 2");
  model2->info = m2info;
  /* Adjust sampling step size */
  abs_sampling_step = args->sampling_step*bbox2_diag;
  abs_sampling_dens = 1/(abs_sampling_step*abs_sampling_step);


  /* Compute the distance from one model to the other */
  dist_surf_surf(model1,model2->mesh,abs_sampling_dens,args->min_sample_freq,
                 &stats,!args->no_gui,(args->quiet ? NULL : progress));

  dist_surf_surf(model2,model1->mesh,abs_sampling_dens,args->min_sample_freq,
                 &stats_rev,0,(args->quiet ? NULL : progress));
  free_face_error(model2->fe);
  model2->fe = NULL;

    /* Print symmetric distance measures */
  output = max(stats.rms_dist,stats_rev.rms_dist)/bbox2_diag*100; 
  return output;
}
