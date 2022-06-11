/**
 * \file sampler-pnt.cpp
 *
 * \brief Samples a PPS built from a PN triangle surface.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date July 2009
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include <exception.hpp>          // DISPLAY_EXCEPTION_INFO()
#include <face_attribute.hpp>     // ppsfrompnt::FaceAttribute
#include <halfedge_attribute.hpp> // ppsfrompnt::HalfedgeAttribute
#include <pps.hpp>                // pps::PPS
#include <ppsfrompnt.hpp>         // ppsfrompnt::PPSfromPNT
#include <ppssampler.hpp>         // pps::PPSsampler
#include <reader.hpp>             // off::Reader
#include <surface.hpp>            // dcel::Surface
#include <vertex_attribute.hpp>   // ppsfrompnt::VertexAttribute
#include <writer.hpp>             // off::Writer

#include <cstdlib>    // atoi, EXIT_SUCCESS, EXIT_FAILURE
#include <filesystem> // std::filesystem
#include <iostream>   // std::cout, std::endl
#include <string>     // std::string

/**
 * \fn int main( int argc , char *argv[ ] )
 *
 * \brief Creates a PPS from a PN triangle. The underlying mesh of the
 * PPS is  given by an input OFF  file, and a PL  approximation to the
 * PPS is written out to another OFF file.
 *
 * \param argc The number of command-line arguments.
 * \param argv The command-line arguments.
 */
int main(int argc, char *argv[]) {
  using std::cout;
  using std::endl;

  namespace fs = std::filesystem;

  using ppsfrompnt::FaceAttribute;
  using ppsfrompnt::HalfedgeAttribute;
  using ppsfrompnt::VertexAttribute;

  using Surface =
      dcel::Surface<VertexAttribute, FaceAttribute, int, HalfedgeAttribute>;

  if (argc != 4) {
    cout << endl;
    cout << "Usage: sampler-pnt <arg 1> <arg 2> <arg 3>" << endl;
    cout << endl;
    cout << "arg1: the full path to the OFF file describing the input triangle "
            "mesh."
         << endl;
    cout << "arg2: the level-of-detail of the output triangular mesh." << endl;
    cout << "arg3: the full path to the output directory." << endl;
    cout << endl;
    cout << "For instance:" << endl;
    cout << "sampler-pnt star.off 3 ./temp" << endl;
    cout << endl;
    
    return EXIT_SUCCESS;
  }

  try {

    /**
     * The  first  command-line  argument   is  the  name  of  the  file
     * describing the  underlying mesh of  the PPS, while the  second is
     * the level of  detail of the PL approximation  to the PPS computed
     * here.
     */
    fs::path inFile(argv[1]);
    unsigned lod = unsigned(atoi(argv[2]));
    fs::path outFolder(argv[3]);

    //
    // Read in the underlying mesh information.
    //

    cout << "Reading input file..." << endl;

    unsigned nv, nf;
    double *vset;
    unsigned *fset;

    off::Reader reader(inFile);

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

    if (vset) {
      delete[] vset;
    }

    if (fset) {
      delete[] fset;
    }

    //
    // Creates a PPS that approximates the PN triangle surface.
    //

    cout << "Creating the PPS..." << endl;

    ppsfrompnt::PPSfromPNT *pps = new ppsfrompnt::PPSfromPNT(mesh);

    pps->build();

    //
    // Sample the PN triangle surface and the PPS.
    //

    cout << "Sampling the PN triangle surface and the PPS ..." << endl;

    pps::PPSsampler<Surface> *sampler = new pps::PPSsampler<Surface>(pps);

    double *lv1;
    double *lv2;

    unsigned *lf;

    sampler->sample(lod, nv, lv1, lv2, nf, lf);

    //
    // Writing out the output files.
    //

    cout << "Generating the output files..." << endl;

    auto filename = inFile.stem();
    const auto suffix = std::to_string(lod) + ".off";
    const auto outFile1 = outFolder / (filename.string() + std::string("-pnt-pps-") + suffix);
    const auto outFile2 = outFolder / (filename.string() + std::string("-pnt-") +  suffix);

    nv /= 3;
    nf /= 3;

    off::Writer writer1(outFile1);
    writer1.write(nv, &lv1[0], nf, &lf[0]);

    if (lv1)
      delete lv1;

    off::Writer writer2(outFile2);
    writer2.write(nv, &lv2[0], nf, &lf[0]);

    if (lv2)
      delete lv2;

    if (lf)
      delete lf;

    //
    // Release memory.
    //

    cout << "Releasing memory..." << endl;

    if (sampler) {
      delete sampler;
    }

    if (pps) {
      delete pps;
    }

    if (mesh) {
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
