#include "fpset.h"

static int enable = 0;

#ifdef MPI

#include <mpi/mpi.h>
#include <stdio.h>

static MPI_File fh;

int fpset_create(const char *fn, int n, int n0)
{
  MPI_Status status;
  int numel, rank, m[2];
  if (MPI_File_open(MPI_COMM_WORLD, fn, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    return -1;
  if (MPI_File_set_size(fh, 0) != MPI_SUCCESS)
    return -1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (!rank) {
    m[0] = n; m[1] = n0;
    if (MPI_File_write_at(fh, 0, m, 2, MPI_INT, &status) != MPI_SUCCESS)
      return -2;
    if (MPI_Get_count(&status, MPI_INT, &numel) != MPI_SUCCESS || numel != 2)
      return -2;
  }
  if (MPI_File_seek_shared(fh, 2 * sizeof(int), MPI_SEEK_SET) != MPI_SUCCESS)
    return -2;
  
  enable = 1;
  return 0;
}

int fpset_open(const char *fn, int n, int n0, int nsamp)
{
  MPI_Status status;
  MPI_Offset fsz, esz;
  int m[2], numel, rank;
    
  if (MPI_File_open(MPI_COMM_WORLD, fn, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    return -1;
  if (MPI_File_get_size(fh, &fsz) != MPI_SUCCESS)
    return -1;

  if (!fsz && !nsamp) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) {
      m[0] = n;
      m[1] = n0;
      if (MPI_File_write_at(fh, 0, m, 2, MPI_INT, &status) != MPI_SUCCESS)
	return -2;
      if (MPI_Get_count(&status, MPI_INT, &numel) != MPI_SUCCESS || numel != 2)
	return -2;
    }
    if (MPI_File_seek_shared(fh, 2 * sizeof(int), MPI_SEEK_SET) != MPI_SUCCESS)
      return -2;
    enable = 1;
    return 0;
  }
  if (MPI_File_read_at_all(fh, 0, m, 2, MPI_INT, &status) != MPI_SUCCESS)
    return -3;
  if (MPI_Get_count(&status, MPI_INT, &numel) != MPI_SUCCESS || numel != 2)
    return -3;

  if (m[0] != n)
    return -4;

  esz = 4 * n * nsamp * sizeof(real_t) + 2 * sizeof(int);
  if (fsz < esz)
    return -5;

  if (MPI_File_seek_shared(fh, esz, MPI_SEEK_SET) != MPI_SUCCESS)
    return -6;

  enable = 1;
  return 0;
}

int fpset_write(int n, real_t *x)
{
  MPI_Status status;
  int numel;
  if (!enable)
    return 0;
  if (MPI_File_write_shared(fh, x, n, MPI_ET_REAL, &status) != MPI_SUCCESS)
    return -1;
  if (MPI_Get_count(&status, MPI_ET_REAL, &numel) != MPI_SUCCESS || numel != n)
    return -1;
  return 0;
}

void fpset_close()
{
  if (enable) {
    MPI_File_close(&fh);
    enable = 0;
  }
}


#else

#include <stdio.h>

static FILE *fh;

int fpset_create(const char *fn, int n, int n0)
{
  int m[2];
  fh = fopen(fn, "wb");
  if (!fh)
    return -1;
  m[0] = n; m[1] = n0;
  if (fwrite(m, sizeof(int), 2, fh) != 2)
    return -2;
  enable = 1;
  return 0;
}

int fpset_open(const char *fn, int n, int n0, int nsamp)
{
  long fsz, esz;
  int m[2];

  fh = fopen(fn, "r+b");
  if (!fh)
    return -1;

  fsz = fread(m, sizeof(int), 2, fh);
  if (!fsz && !nsamp) {
    m[0] = n; m[1] = n0;
    if (fwrite(m, sizeof(int), 2, fh) != 2)
      return -2;
    enable = 1;
    return 0;
  }
  if (fsz != 2)
    return -3;

  if (m[0] != n)
    return -4;

  if (!fseek(fh, 0, SEEK_END))
    return -6;
  fsz = ftell(fh);
  if (fsz < 0)
    return -6;

  esz = 4 * n * nsamp * sizeof(real_t) + 2 * sizeof(int);
  if (fsz < esz)
    return -5;

  if (!fseek(fh, esz, SEEK_SET))
    return -6;

  enable = 1;
  return 0;
}

int fpset_write(int n, real_t *x)
{
  if (!enable)
    return 0;
  return fwrite(x, sizeof(real_t), n, fh) != n;
}

void fpset_close()
{
  if (enable) {
    fclose(fh);
    enable = 0;
  }
}


#endif
