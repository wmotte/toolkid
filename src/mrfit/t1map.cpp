#include "mrfitCommon.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "procparser.h"

int main(int argc, char ** argv) {

	tkd::CmdParser p("t1map", "T1 mapping");

	std::string inputFileName;
	std::string outputFileName;
	std::string procparFileName;
	std::string outputFileNameA;
	std::string outputFileNameC;
	std::string outputFileNameR1;

	int algorithm = 0;
	double maxT1 = 10000.0;
	std::vector<double> times;

	p.AddArgument(inputFileName, "input")
		->AddAlias("i")
		->SetDescription("Input 4D image (3D x number of times)")
		->SetRequired(true);

	p.AddArgument(procparFileName, "procpar")
		->AddAlias("p")
		->SetDescription("Path to FID/procpar");

	p.AddArgument(times, "times")
		->AddAlias("t")
		->SetDescription("TR times (in seconds)");

	p.AddArgument(outputFileName, "output")
		->AddAlias("o")
		->SetDescription("Output T1 map")
		->SetRequired(true);

	p.AddArgument(outputFileNameA, "output-A")
		->AddAlias("oA")
		->SetDescription("Output map of fitted constant A");

	p.AddArgument(outputFileNameC, "output-C")
		->AddAlias("oC")
		->SetDescription("Output map of fitted constant C");

	p.AddArgument(outputFileNameR1, "output-R1")
		->AddAlias("oR1")
		->SetDescription("Output map of R^2 fitting values");

	p.AddArgument(algorithm, "algorithm")
		->AddAlias("a")
		->SetDescription("Fitting algorithm: "
				"0=IDEAL_STEADY_STATE, "
				"1=INVERSION_RECOVERY, "
				"2=ABSOLUTE_INVERSION_RECOVERY, "
				"3=LOOK_LOCKER, "
				"4=ABSOLUTE_LOOK_LOCKER, "
				"5=HYBRID_STEADY_STATE_3PARAM, "
				"6=INVERSION_RECOVERY_3PARAM, "
				"7=ABSOLUTE_INVERSION_RECOVERY_3PARAM, "
				"8=ABSOLUTE_LOOK_LOCKER_3PARAM"
				"(default: 3)");


	p.AddArgument(maxT1, "maximum-T1")
		->AddAlias("max")
		->SetDescription("Maximum T1 (default: 10000)");

	if (!p.Parse(argc, argv)) {
		p.PrintUsage(std::cout);
		return -1;
	}

	Procparser pp;

	if (times.size() == 0) {
		if (procparFileName == "") {
			std::cerr << "Either provide the path to the FID's procpar, or supply the times" << std::endl;
			return -1;
		}

		pp.Parse(procparFileName);
		if (pp.Has("TR")) {
			for (int i = 0; i < pp.GetSize("TR"); ++i) {
				times.push_back(pp.GetAs<double> ("TR", i) / 1000.);
			}
		} else { //TODO
			for (int i = 1; i <= pp.GetAs<int> ("ne"); ++i) {
				times.push_back(pp.GetAs<double> ("tr") * static_cast<double> (i));
			}
		}
	}

	typedef itk::Image<float, 4> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inputFileName.c_str());
	reader->Update();

	mrfit::MRFit::Pointer fit = mrfit::MRFit::New();
	fit->SetInput(reader->GetOutput());
	fit->SetTimes(times);
	fit->FitT1(static_cast<mrfit::MRFit::T1FittingType> (algorithm), maxT1);

	typedef itk::Image<float, 3> OutputImageType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;

	WriterType::Pointer writer = WriterType::New();

	// output T1
	writer->SetInput(fit->GetMap(0));
	writer->SetFileName(outputFileName.c_str());
	writer->Update();

	// output A
	if (outputFileNameA != "") {
		writer->SetInput(fit->GetMap(1));
		writer->SetFileName(outputFileNameA.c_str());
		writer->Update();
	}

	// output C
	if (outputFileNameC != "") {
		writer->SetInput(fit->GetMap(2));
		writer->SetFileName(outputFileNameC.c_str());
		writer->Update();
	}

	// output r^2
	if (outputFileNameR1 != "") {
		writer->SetInput(fit->GetMap(3));
		writer->SetFileName(outputFileNameR1.c_str());
		writer->Update();
	}

	return 0;
}
