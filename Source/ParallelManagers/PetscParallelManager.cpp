// RealType*                 bufferleft[sizeY];
// PressureBufferFillStencil pfill(params, buffers)
// PressureBufferReadStencil pread(params, buffers)
// similar for vels

// communicate();
// pfill.iterate();
// velsfill.iterate(); for
//   u and v

//     MPI_sendrecv();

// pread.iterate();
// velread.iterate();