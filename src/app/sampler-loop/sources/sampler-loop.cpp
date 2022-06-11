/**
 * \file sampler-loop.cpp
 *
 * \brief Samples a PPS built from a Loop subdivision surface.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date January 2010
 */

#include <exception.hpp>          // DISPLAY_EXCEPTION_INFO()
#include <face_attribute.hpp>     // ppsfromloop::FaceAttribute
#include <halfedge_attribute.hpp> // ppsfromloop::HalfedgeAttribute
#include <pps.hpp>                // pps::PPS
#include <ppsfromloop.hpp>        // ppsfromloop::PPSfromLOOP
#include <ppssampler.hpp>         // pps::PPSsampler
#include <reader.hpp>             // off::Reader
#include <surface.hpp>            // dcel::Surface
#include <vertex_attribute.hpp>   // ppsfromloop::VertexAttribute
#include <writer.hpp>             // off::Writer

#include <cstdlib>    // atoi, EXIT_SUCCESS, EXIT_FAILURE
#include <filesystem> // std::filesystem
#include <iostream>   // std::cout, std::endl
#include <string>     // std::string

/**
 * \fn int main( int argc , char *argv[ ] )
 *
 * \brief  Creates  a  PPS  from  a  Loop  subdivision  surface.   The
 * underlying mesh  of the PPS  is given by an  input OFF file,  and a
 * piecewise-linear approximation to the PPS is written out to another
 * OFF file.
 *
 * \param argc The number of command-line arguments.
 * \param argv The command-line arguments.
 */
int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;

  namespace fs = std::filesystem;

  using ppsfromloop::FaceAttribute;
  using ppsfromloop::HalfedgeAttribute;
  using ppsfromloop::VertexAttribute;

  using Surface =
      dcel::Surface<VertexAttribute, FaceAttribute, int, HalfedgeAttribute>;

  if (argc != 5) {
    cout << endl;
    cout << "Usage: sampler-loop <arg 1> <arg 2> <arg 3> <arg 4>" << endl;
    cout << endl;
    cout << "arg1: the full path to the OFF file describing the input triangle "
            "mesh."
         << endl;
    cout << "arg2: the level-of-detail of the output triangular mesh." << endl;
    cout << "arg3: the name of the Loop surface evaluator data table file."
         << endl;
    cout << "arg4: the full path to the output directory." << endl;
    cout << endl;
    cout << "For instance:" << endl;
    cout << "sampler-loop star 3 lpdata50NT.dat ./temp" << endl;
    cout << endl;

    return EXIT_SUCCESS;
  }

  try {

    /**
     * The  first  command-line  argument   is  the  name  of  the  file
     * describing the underlying mesh of the PPS, the second argument is
     * the level of  detail of the PL approximation  to the PPS computed
     * here, and the  third argument is the name  of the file containing
     * the Loop surface evaluator data table.
     */
    fs::path inMeshFile(argv[1]);
    unsigned lod = unsigned(atoi(argv[2]));
    std::string inEvalFile(argv[3]);
    fs::path outFolder(argv[4]);

    //
    // Read in the underlying mesh information.
    //

    cout << "Reading input file..." << endl;

    unsigned nv, nf;
    double *vset;
    unsigned *fset;

    off::Reader reader(inMeshFile);

    reader.read(nv, vset, nf, fset);

    //
    // Creates the underlyng mesh and stores it in a DCEL.
    //

    cout << "Creating surface..." << endl;

    Surface *mesh = new Surface(nv, vset, nf, fset);

    //
    // Release memory.
    //

    cout << "Releasing memory..." << endl;

    if (vset != 0) {
      delete[] vset;
    }

    if (fset != 0) {
      delete[] fset;
    }

    //
    // Creates a PPS that approximates the Loop subdivision surface.
    //

    cout << "Creating the PPS..." << endl;

    ppsfromloop::PPSfromLOOP *pps =
        new ppsfromloop::PPSfromLOOP(mesh, inEvalFile);

    pps->build();

    //
    // Sample the Loop subdivision surface and the PPS.
    //

    cout << "Sampling the Loop subdivision surface and the PPS ..." << endl;

    pps::PPSsampler<Surface> *sampler = new pps::PPSsampler<Surface>(pps);

    double *lv1;
    double *lv2;

    unsigned *lf;

    sampler->sample(lod, nv, lv1, lv2, nf, lf);

    //
    // Writing out the output files.
    //

    cout << "Generating the output files..." << endl;

    auto filename = inMeshFile.stem();
    const auto suffix = std::to_string(lod) + ".off";
    const auto outFile1 =
        outFolder / (filename.string() + std::string("-loop-pps-") + suffix);
    const auto outFile2 =
        outFolder / (filename.string() + std::string("-loop-") + suffix);

    nv /= 3;
    nf /= 3;

    off::Writer writer1(outFile1);
    writer1.write(nv, &lv1[0], nf, &lf[0]);

    if (lv1 != 0)
      delete lv1;

    off::Writer writer2(outFile2);
    writer2.write(nv, &lv2[0], nf, &lf[0]);

    if (lv2 != 0)
      delete lv2;

    if (lf != 0)
      delete lf;

    //
    // Release memory.
    //

    cout << "Releasing memory..." << endl;

    if (sampler != 0) {
      delete sampler;
    }

    if (pps != 0) {
      delete pps;
    }

    if (mesh != 0) {
      delete mesh;
    }

  } catch (const utils::Exception &xpt) {
    std::cout << std::endl;
    DISPLAY_EXCEPTION_INFO(xpt);
    std::cout << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
