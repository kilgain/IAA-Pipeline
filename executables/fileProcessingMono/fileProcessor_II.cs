// --------------------------------------------------------------------------------------------------------------------
// <copyright file="Program.cs" company="Thermo Fisher Scientific">
//   Copyright © Thermo Fisher Scientific. All rights reserved.
// </copyright>
// This source file contains the code for the console application that demonstrates who to
// read our RAW files using the RAWFileReader .Net library.  This example uses only a fraction
// read our RAW files using the RAWFileReader .Net library.  This example uses only a fraction
// of the functions available in the RawFileReader library.  Please consult the RawFileReader
// documentation for a complete list of methods available.
// 
// This code has been compiled and tested using Visual Studio 2013 Update 5, Visual Studio 
// 2015 Update 3, Visual Studio 2017 Update 2, Visual Studio MAC, and MonoDevelop.  Additional 
// tools used include Resharper 2017.1.  This application has been tested with .Net Framework 4.5.1, 
// 4.5.2, 4.6, 4.6.1, 4.6.2, and 4.7.
//
// Questions about the program can be directed to jim.shofstahl@thermofisher.com
// --------------------------------------------------------------------s------------------------------------------------
namespace RawFileReaderDotNetExample
{
    // These are the libraries necessary to read the Thermo RAW files.  The Interfaces
    // library contains the extension for accessing the scan averaging and background
    // subtraction methods.
    using System;
    using System.Collections.Generic;
    using System.Diagnostics;
    using System.IO;
    using System.Runtime.ExceptionServices;

    using ThermoFisher.CommonCore.Data;
    using ThermoFisher.CommonCore.Data.Business;
    using ThermoFisher.CommonCore.Data.FilterEnums;
    using ThermoFisher.CommonCore.Data.Interfaces;
    using ThermoFisher.CommonCore.MassPrecisionEstimator;
    using ThermoFisher.CommonCore.RawFileReader;


    /// <summary>
    /// A C# example program showing how to use RAWFileReader.  More information on the RAWFileReader methods used
    /// in this example and the other methods available in RAWFileReader can be found in the RAWFileReader user
    /// documentation, that is installed with the RAWFileReader software.
    /// This program has been tested with RAWFileReader 4.0.22  Changes maybe necessary with other versions
    /// of RAWFileReader.
    /// </summary>
    internal static class Program
    {

        /// <summary>
        /// The main routine for this example program.  This example program only shows how to use the RawFileReader library
        /// in a single-threaded application but our documentation for RawFileReader describes how to use it in a multi-threaded
        /// application.
        /// </summary>
        /// <param name="args">The command line arguments for this program.  The RAW file name should be passed as the first argument</param>
        private static void Main(string[] args)
        {
            // This local variable controls if the AnalyzeAllScans method is called
            bool analyzeScans = false;


            ///allen add
            //FileStream ostrm;
            //StreamWriter writer;
            //TextWriter oldOut = Console.Out;
            //try
            //{
            //    ostrm = new FileStream("./Redirect.txt", FileMode.OpenOrCreate, FileAccess.Write);
            //    writer = new StreamWriter(ostrm);
            //}
            //catch (Exception e)
            //{
            //    Console.WriteLine("Cannot open Redirect.txt for writing");
            //    Console.WriteLine(e.Message);
             //   return;
            //}
            //Console.SetOut(writer);
            //Console.WriteLine("This is a line of text");
            //Console.WriteLine("Everything written to Console.Write() or");
            //Console.WriteLine("Console.WriteLine() will be written to a file");
     
            //writer.Close();
            //ostrm.Close();
            //Console.WriteLine("Done");

            ///allen stop 
            ///
            // Get the memory used at the beginning of processing
            Process processBefore = Process.GetCurrentProcess();
            long memoryBefore = processBefore.PrivateMemorySize64 / 1024;

            try
            { 
                // Check to see if the RAW file name was supplied as an argument to the program
                string filename = string.Empty;

                if (args.Length > 0)
                {
                    filename = args[0];
                }

                if (string.IsNullOrEmpty(filename))
                {
                    Console.WriteLine("No RAW file specified!");

                    return;
                }


                    FileInfo fi = new FileInfo(filename);
                    string Basename = fi.Name;


                    //System.Console.Write("{0} ", i);
                    //System.Console.Write("is a file");

                    // Create the IRawDataPlus object for accessing the RAW file
                    var rawFile = RawFileReaderAdapter.FileFactory(filename);
			
                    Console.WriteLine(filename);
                    Console.WriteLine(rawFile);
                    if (!rawFile.IsOpen || rawFile.IsError)
                    {
                        Console.WriteLine("Unable to access the RAW file using the RawFileReader class!");

                        return;
                    }

                    // Check for any errors in the RAW file
                    if (rawFile.IsError)
                    {
                        Console.WriteLine("Error opening ({0}) - {1}", rawFile.FileError, filename);

                        return;
                    }

                    // Check if the RAW file is being acquired
                    if (rawFile.InAcquisition)
                    {
                        Console.WriteLine("RAW file still being acquired - " + filename);

                        return;
                    }

                    // Get the number of instruments (controllers) present in the RAW file and set the 
                    // selected instrument to the MS instrument, first instance of it
                    Console.WriteLine("The RAW file has data from {0} instruments" + rawFile.InstrumentCount);

                    rawFile.SelectInstrument(Device.MS, 1);

                    // Get the first and last scan from the RAW file
                    int firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum;
                    int lastScanNumber = rawFile.RunHeaderEx.LastSpectrum;

                    // Get the start and end time from the RAW file
                    double startTime = rawFile.RunHeaderEx.StartTime;
                    double endTime = rawFile.RunHeaderEx.EndTime;

                    // Print some OS and other information
                    Console.WriteLine("System Information:");
                    Console.WriteLine("   OS Version: " + Environment.OSVersion);
                    Console.WriteLine("   64 bit OS: " + Environment.Is64BitOperatingSystem);
                    Console.WriteLine("   Computer: " + Environment.MachineName);
                    Console.WriteLine("   # Cores: " + Environment.ProcessorCount);
                    Console.WriteLine("   Date: " + DateTime.Now);
                    Console.WriteLine();

                    // Get some information from the header portions of the RAW file and display that information.
                    // The information is general information pertaining to the RAW file.
                    Console.WriteLine("General File Information:");
                    Console.WriteLine("   RAW file: " + rawFile.FileName);
                    Console.WriteLine("   RAW file version: " + rawFile.FileHeader.Revision);
                    Console.WriteLine("   Creation date: " + rawFile.FileHeader.CreationDate);
                    Console.WriteLine("   Operator: " + rawFile.FileHeader.WhoCreatedId);
                    Console.WriteLine("   Number of instruments: " + rawFile.InstrumentCount);
                    Console.WriteLine("   Description: " + rawFile.FileHeader.FileDescription);
                    Console.WriteLine("   Instrument model: " + rawFile.GetInstrumentData().Model);
                    Console.WriteLine("   Instrument name: " + rawFile.GetInstrumentData().Name);
                    Console.WriteLine("   Serial number: " + rawFile.GetInstrumentData().SerialNumber);
                    Console.WriteLine("   Software version: " + rawFile.GetInstrumentData().SoftwareVersion);
                    Console.WriteLine("   Firmware version: " + rawFile.GetInstrumentData().HardwareVersion);
                    Console.WriteLine("   Units: " + rawFile.GetInstrumentData().Units);
                    Console.WriteLine("   Mass resolution: {0:F3} ", rawFile.RunHeaderEx.MassResolution);
                    Console.WriteLine("   Number of scans: {0}", rawFile.RunHeaderEx.SpectraCount);
                    Console.WriteLine("   Scan range: {0} - {1}", firstScanNumber, lastScanNumber);
                    Console.WriteLine("   Time range: {0:F2} - {1:F2}", startTime, endTime);
                    Console.WriteLine("   Mass range: {0:F4} - {1:F4}", rawFile.RunHeaderEx.LowMass, rawFile.RunHeaderEx.HighMass);
                    Console.WriteLine();

                    // Get information related to the sample that was processed
                    Console.WriteLine("Sample Information:");
                    Console.WriteLine("   Sample name: " + rawFile.SampleInformation.SampleName);
                    Console.WriteLine("   Sample id: " + rawFile.SampleInformation.SampleId);
                    Console.WriteLine("   Sample type: " + rawFile.SampleInformation.SampleType);
                    Console.WriteLine("   Sample comment: " + rawFile.SampleInformation.Comment);
                    Console.WriteLine("   Sample vial: " + rawFile.SampleInformation.Vial);
                    Console.WriteLine("   Sample volume: " + rawFile.SampleInformation.SampleVolume);
                    Console.WriteLine("   Sample injection volume: " + rawFile.SampleInformation.InjectionVolume);
                    Console.WriteLine("   Sample row number: " + rawFile.SampleInformation.RowNumber);
                    Console.WriteLine("   Sample dilution factor: " + rawFile.SampleInformation.DilutionFactor);
                    Console.WriteLine();

                    // Read the first instrument method (most likely for the MS portion of the instrument).
                    // NOTE: This method reads the instrument methods from the RAW file but the underlying code
                    // uses some Microsoft code that hasn't been ported to Linux or MacOS.  Therefore this
                    // method won't work on those platforms therefore the check for Windows.
                    if (Environment.OSVersion.ToString().Contains("Windows"))
                    {
                        var deviceNames = rawFile.GetAllInstrumentNamesFromInstrumentMethod();

                        foreach (var device in deviceNames)
                        {
                            Console.WriteLine("Instrument method: " + device);
                        }

                        Console.WriteLine();
                    }

                    // Display all of the trailer extra data fields present in the RAW file
                    ListTrailerExtraFields(rawFile);

                    // Get the number of filters present in the RAW file
                    int numberFilters = rawFile.GetFilters().Count;

                    // Get the scan filter for the first and last spectrum in the RAW file
                    var firstFilter = rawFile.GetFilterForScanNumber(firstScanNumber);
                    var lastFilter = rawFile.GetFilterForScanNumber(lastScanNumber);

                    Console.WriteLine("Filter Information:");
                    Console.WriteLine("   Scan filter (first scan): " + firstFilter.ToString());
                    Console.WriteLine("   Scan filter (last scan): " + lastFilter.ToString());
                    Console.WriteLine("   Total number of filters:" + numberFilters);
                    Console.WriteLine();

                    // Get the BasePeak chromatogram for the MS data
                    GetChromatogram(rawFile, firstScanNumber, lastScanNumber, true);




                    var myoutname = args[1] + "/" + Basename + "output.txt";


                    // Read the scan information for each scan in the RAW file
                    ReadScanInformation(rawFile, firstScanNumber, lastScanNumber, true, myoutname);

                    //for the hell of it, set this to false
                    //ReadScanInformation(rawFile, firstScanNumber, lastScanNumber, false);

                    // Get a spectrum from the RAW file.  
                    GetSpectrum(rawFile, firstScanNumber, firstFilter.ToString(), false);

                    // Get a spectrum from the RAW file.  , allen add: set boolean to true for the hell of it
                    GetSpectrum(rawFile, firstScanNumber, firstFilter.ToString(), true);

                    // Get a average spectrum from the RAW file for the first 15 scans in the file.  
                    GetAverageSpectrum(rawFile, 1, 15, false);


                    //rerun with true, allen edit - see what happens
                    GetAverageSpectrum(rawFile, 1, 15, true);


                    // Read each spectrum
                    ReadAllSpectra(rawFile, firstScanNumber, lastScanNumber, true);

                    // Calculate the mass precision for a spectrum
                    CalculateMassPrecision(rawFile, 1);


                    // Check all of the scans for out of order data.  This method isn't enabled by
                    // default because it is very, very time consuming.  If you would like to 
                    // call this method change the value of _analyzeScans to true.
                    if (analyzeScans)
                    {
                        AnalyzeAllScans(rawFile, firstScanNumber, lastScanNumber);
                    }

                    // Close (dispose) the RAW file
                    Console.WriteLine();
                    Console.WriteLine("Closing " + filename);

                    rawFile.Dispose();
                //}
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error accessing RAWFileReader library! - " + ex.Message);
            }



//            Console.SetOut(oldOut);
//            writer.Close();
//            ostrm.Close();
//            Console.WriteLine("Done");
            // Get the memory used at the end of processing
            Process processAfter = Process.GetCurrentProcess();
            long memoryAfter = processAfter.PrivateMemorySize64 / 1024;

            Console.WriteLine();
            Console.WriteLine("Memory Usage:");
            Console.WriteLine("   Before {0} kb, After {1} kb, Extra {2} kb", memoryBefore, memoryAfter, memoryAfter - memoryBefore);


            string path = Directory.GetCurrentDirectory();
            //string target = @"c:\temp";
            Console.WriteLine("The current directory is {0}", path);

          

        }
        
        /// <summary>
        /// Reads all of the scans in the RAW and looks for out of order data
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan in the RAW file
        /// </param>
        /// <param name="lastScanNumber">
        /// the last scan in the RAW file
        /// </param>
        private static void AnalyzeAllScans(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber)
        {

            Console.WriteLine("This is a line of text");
            Console.WriteLine("Everything written to Console.Write() or");
            Console.WriteLine("Console.WriteLine() will be written to a file");

            


            // Test the preferred (normal) data and centroid (high resolution/label) data
            int failedCentroid = 0;
            int failedPreferred = 0;

            for (int scanNumber = firstScanNumber; scanNumber <= lastScanNumber; scanNumber++)
            {
                // Get each scan from the RAW file
                var scan = Scan.FromFile(rawFile, scanNumber);

                // Check to see if the RAW file contains label (high-res) data and if it is present
                // then look for any data that is out of order
                if (scan.HasCentroidStream)
                {
                    if (scan.CentroidScan.Length > 0)
                    {
                        double currentMass = scan.CentroidScan.Masses[0];

                        for (int index = 1; index < scan.CentroidScan.Length; index++)
                        {
                            if (scan.CentroidScan.Masses[index] > currentMass)
                            {
                                currentMass = scan.CentroidScan.Masses[index];
                            }
                            else
                            {
                                if (failedCentroid == 0)
                                {
                                    Console.WriteLine("First failure: Failed in scan data at: Scan: " + scanNumber + " Mass: "
                                        + currentMass.ToString("F4"));
                                }

                                failedCentroid++;
                            }
                        }
                    }
                }

                // Check the normal (non-label) data in the RAW file for any out-of-order data
                if (scan.PreferredMasses.Length > 0)
                {
                    double currentMass = scan.PreferredMasses[0];

                    for (int index = 1; index < scan.PreferredMasses.Length; index++)
                    {
                        if (scan.PreferredMasses[index] > currentMass)
                        {
                            currentMass = scan.PreferredMasses[index];
                        }
                        else
                        {
                            if (failedPreferred == 0)
                            {
                                Console.WriteLine("First failure: Failed in scan data at: Scan: " + scanNumber + " Mass: "
                                    + currentMass.ToString("F2"));
                            }

                            failedPreferred++;
                        }
                    }
                }
            }

            // Display a message indicating if any of the scans had data that was "out of order"
            if (failedPreferred == 0 && failedCentroid == 0)
            {
                Console.WriteLine("Analysis completed: No out of order data found");
            }
            else
            {
                Console.WriteLine("Analysis completed: Preferred data failed: " + failedPreferred + " Centroid data failed: " + failedCentroid);
            }

            //Console.SetOut(oldOut);
            //writer.Close();
            //ostrm.Close();
            Console.WriteLine("Done");

        }

        /// <summary>
        /// Calculates the mass precision for a spectrum.  
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="scanNumber">
        /// The scan to process
        /// </param>
        private static void CalculateMassPrecision(IRawDataPlus rawFile, int scanNumber)
        {
            // Get the scan from the RAW file
            var scan = Scan.FromFile(rawFile, scanNumber);

            // Get the scan event and from the scan event get the analyzer type for this scan
            var scanEvent = rawFile.GetScanEventForScanNumber(scanNumber);

            // Get the trailer extra data to get the ion time for this file
            LogEntry logEntry = rawFile.GetTrailerExtraInformation(scanNumber);

            var trailerHeadings = new List<string>();
            var trailerValues = new List<string>();
            for (var i = 0; i < logEntry.Length; i++)
            {
                trailerHeadings.Add(logEntry.Labels[i]);
                trailerValues.Add(logEntry.Values[i]);
            }

            // Create the mass precision estimate object
            IPrecisionEstimate precisionEstimate = new PrecisionEstimate();

            // Get the ion time from the trailer extra data values
            var ionTime = precisionEstimate.GetIonTime(scanEvent.MassAnalyzer, scan, trailerHeadings, trailerValues);

            // Calculate the mass precision for the scan
            var listResults = precisionEstimate.GetMassPrecisionEstimate(scan, scanEvent.MassAnalyzer, ionTime, rawFile.RunHeader.MassResolution);

            // Output the mass precision results
            if (listResults.Count > 0)
            {
                Console.WriteLine("Mass Precision Results:");

                foreach (var result in listResults)
                {
                    Console.WriteLine("Mass {0:F5}, mmu = {1:F3}, ppm = {2:F2}", result.Mass, result.MassAccuracyInMmu, result.MassAccuracyInPpm);
                }
            }
        }

        /// <summary>
        /// Gets the average spectrum from the RAW file.  
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan to consider for the averaged spectrum
        /// </param>
        /// <param name="lastScanNumber">
        /// The last scan to consider for the averaged spectrum
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetAverageSpectrum(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData)
        {

            ///allen add
            //FileStream ostrm;
            //StreamWriter writer;
            //TextWriter oldOut = Console.Out;
            //try
            //{
            //    ostrm = new FileStream("./RedirectIII.txt", FileMode.OpenOrCreate, FileAccess.Write);
            //    writer = new StreamWriter(ostrm);
            //}
            //catch (Exception e)
            //{
            //    Console.WriteLine("Cannot open RedirectIII.txt for writing");
            //    Console.WriteLine(e.Message);
            //    return;
            //}
            //Console.SetOut(writer);
            //Console.WriteLine("This is a line of text");
            //Console.WriteLine("Everything written to Console.Write() or");
            //Console.WriteLine("Console.WriteLine() will be written to a file");


            // Create the mass options object that will be used when averaging the scans
            var options = rawFile.DefaultMassOptions();

            options.ToleranceUnits = ToleranceUnits.ppm;
            options.Tolerance = 5.0;

            // Get the scan filter for the first scan.  This scan filter will be used to located
            // scans within the given scan range of the same type
            var scanFilter = rawFile.GetFilterForScanNumber(firstScanNumber);

            // Get the average mass spectrum for the provided scan range. In addition to getting the
            // average scan using a scan range, the library also provides a similar method that takes
            // a time range.
            var averageScan = rawFile.AverageScansInScanRange(firstScanNumber, lastScanNumber, scanFilter, options);

            if (averageScan.HasCentroidStream)
            {
                Console.WriteLine("Average spectrum ({0} points)", averageScan.CentroidScan.Length);

                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < averageScan.CentroidScan.Length; i++)
                    {
                       Console.WriteLine("what we want");
                       Console.WriteLine("  {0:F4} {1:F0}", averageScan.CentroidScan.Masses[0], averageScan.CentroidScan.Intensities[i]);
                    }
                }
            }

            // This example uses a different method to get the same average spectrum that was calculated in the
            // previous portion of this method.  Instead of passing the start and end scan, a list of scans will
            // be passed to the GetAveragedMassSpectrum function.
            List<int> scans = new List<int>(new[] { 1, 6, 7, 9, 11, 12, 14 });

            averageScan = rawFile.AverageScans(scans, options);

            if (averageScan.HasCentroidStream)
            {
                Console.WriteLine("Average spectrum ({0} points)", averageScan.CentroidScan.Length);
                Console.WriteLine("going into centroid stream");


                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < averageScan.CentroidScan.Length; i++)
                    {
                        Console.WriteLine("  {0:F4} {1:F0}", averageScan.CentroidScan.Masses[0], averageScan.CentroidScan.Intensities[i]);
                        Console.WriteLine("within the centroid stream");
                    }
                }
            }

            //Console.WriteLine();

            //Console.SetOut(oldOut);
            //writer.Close();
            //ostrm.Close();
            //Console.WriteLine("Done");

        }

        /// <summary>
        /// Reads the base peak chromatogram for the RAW file
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="startScan">
        /// Start scan for the chromatogram
        /// </param>
        /// <param name="endScan">
        /// End scan for the chromatogram
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetChromatogram(IRawDataPlus rawFile, int startScan, int endScan, bool outputData)
        {
            // Define the settings for getting the Base Peak chromatogram
            ChromatogramTraceSettings settings = new ChromatogramTraceSettings(TraceType.BasePeak);

            // Get the chromatogram from the RAW file. 
            var data = rawFile.GetChromatogramData(new IChromatogramSettings[] { settings }, startScan, endScan);

            // Split the data into the chromatograms
            var trace = ChromatogramSignal.FromChromatogramData(data);

            if (trace[0].Length > 0)
            {
                // Print the chromatogram data (time, intensity values)
                Console.WriteLine("Base Peak chromatogram ({0} points)", trace[0].Length);

                if (outputData)
                {
                    for (int i = 0; i < trace[0].Length; i++)
                    {
                        Console.WriteLine("  {0} - {1:F3}, {2:F0}", i, trace[0].Times[i], trace[0].Intensities[i], trace[0].SignalScans);
                        Console.WriteLine("  {0} - {1:F3}, {2:F0}", i, trace[0].Times[i], trace[0].Intensities[i], trace[0].SignalScans);
                        Console.WriteLine("more basepeak stuff");
                    }
                }
            }

            //Console.WriteLine();
        }

        /// <summary>
        /// Gets the spectrum from the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="scanNumber">
        /// The scan number being read
        /// </param>
        /// <param name="scanFilter">
        /// The scan filter for that scan
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void GetSpectrum(IRawDataPlus rawFile, int scanNumber, string scanFilter, bool outputData)
        {
            // Check for a valid scan filter
            if (string.IsNullOrEmpty(scanFilter))
            {
                return;
            }

            // Get the scan statistics from the RAW file for this scan number
            var scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber);

            // Check to see if the scan has centroid data or profile data.  Depending upon the
            // type of data, different methods will be used to read the data.  While the ReadAllSpectra
            // method demonstrates reading the data using the Scan.FromFile method, generating the
            // Scan object takes more time and memory to do, so that method isn't optimum.
            if (scanStatistics.IsCentroidScan)
            {
                // Get the centroid (label) data from the RAW file for this scan
                var centroidStream = rawFile.GetCentroidStream(scanNumber, false);

                Console.WriteLine("Spectrum (centroid/label) {0} - {1} points", scanNumber, centroidStream.Length);

                // Print the spectral data (mass, intensity, charge values).  Not all of the information in the high resolution centroid 
                // (label data) object is reported in this example.  Please check the documentation for more information about what is
                // available in high resolution centroid (label) data.
                if (outputData)
                {
                    for (int i = 0; i < centroidStream.Length; i++)
                    {
                        Console.WriteLine("  {0} - {1:F4}, {2:F0}, {3:F0}", i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]);
                    }
                }

                Console.WriteLine();
            }
            else
            {
                // Get the segmented (low res and profile) scan data
                var segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics);

                Console.WriteLine("Spectrum (normal data) {0} - {1} points", scanNumber, segmentedScan.Positions.Length);

                // Print the spectral data (mass, intensity values)
                if (outputData)
                {
                    for (int i = 0; i < segmentedScan.Positions.Length; i++)
                    {
                        Console.WriteLine("  {0} - {1:F4}, {2:F0}", i, segmentedScan.Positions[i], segmentedScan.Intensities[i]);

                    }
                    Console.WriteLine("profile data");
                }

                Console.WriteLine();
            }
        }

        /// <summary>
        /// Reads and reports the trailer extra data fields present in the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file
        /// </param>
        private static void ListTrailerExtraFields(IRawDataPlus rawFile)
        {
            // Get the Trailer Extra data fields present in the RAW file
            var trailerFields = rawFile.GetTrailerExtraHeaderInformation();

            // Display each value
            int i = 0;

            Console.WriteLine("Trailer Extra Data Information:");

            foreach (var field in trailerFields)
            {
                Console.WriteLine("   Field {0} = {1} storing data of type {2}", i, field.Label, field.DataType);

                i++;
            }

            Console.WriteLine();
        }

        /// <summary>
        /// Read all spectra in the RAW file.
        /// </summary>
        /// <param name="rawFile">
        /// The raw file.
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan number.
        /// </param>
        /// <param name="lastScanNumber">
        /// The last scan number.
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        [HandleProcessCorruptedStateExceptions]
        private static void ReadAllSpectra(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData)
        {


            ///allen add
            //FileStream ostrm;
            //StreamWriter writer;
            //TextWriter oldOut = Console.Out;
            //try
            //{
            //    ostrm = new FileStream("./RedirectIV.txt", FileMode.OpenOrCreate, FileAccess.Write);
            //    writer = new StreamWriter(ostrm);
            //}
            //catch (Exception e)
            //{
             //   Console.WriteLine("Cannot open RedirectIV.txt for writing");
            //    Console.WriteLine(e.Message);
            //    return;
            //}
            //Console.SetOut(writer);

            for (int scanNumber = firstScanNumber; scanNumber <= lastScanNumber; scanNumber++)
            {
                try
                {
                    // Get the scan filter for the spectrum
                    var scanFilter = rawFile.GetFilterForScanNumber(firstScanNumber);  
                    
                    if (string.IsNullOrEmpty(scanFilter.ToString()))
                    {
                        continue;
                    }

                    // Get the scan from the RAW file.  This method uses the Scan.FromFile method which returns a
                    // Scan object that contains both the segmented and centroid (label) data from an FTMS scan
                    // or just the segmented data in non-FTMS scans.  The GetSpectrum method demonstrates an
                    // alternative method for reading scans.
                    var scan = Scan.FromFile(rawFile, scanNumber);
                    
                    // If that scan contains FTMS data then Centroid stream will be populated so check to see if it is present.
                    int labelSize = 0;

                    if (scan.HasCentroidStream)
                    {
                        labelSize = scan.CentroidScan.Length;
                    }

                    // For non-FTMS data, the preferred data will be populated
                    int dataSize = scan.PreferredMasses.Length;

                    if (outputData)
                    {
                        Console.WriteLine("Spectrum {0} - {1}: normal {2}, label {3} points", scanNumber, scanFilter.ToString(), dataSize, labelSize);
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Error reading spectrum {0} - {1}", scanNumber, ex.Message);
                }
            }

            //Console.SetOut(oldOut);
            //writer.Close();
            //ostrm.Close();
            //Console.WriteLine("Done");

        }

        /// <summary>
        /// Reads the general scan information for each scan in the RAW file using the scan filter object and also the
        /// trailer extra data section for that same scan.
        /// </summary>
        /// <param name="rawFile">
        /// The RAW file being read
        /// </param>
        /// <param name="firstScanNumber">
        /// The first scan in the RAW file
        /// </param>
        /// <param name="lastScanNumber">
        /// the last scan in the RAW file
        /// </param>
        /// <param name="outputData">
        /// The output data flag.
        /// </param>
        private static void ReadScanInformation(IRawDataPlus rawFile, int firstScanNumber, int lastScanNumber, bool outputData, string myoutname)
        {
            
            
            
            ///allen add
            FileStream ostrm;
            StreamWriter writer;
            TextWriter oldOut = Console.Out;
            try
            {
                ostrm = new FileStream(myoutname, FileMode.OpenOrCreate, FileAccess.Write);
                writer = new StreamWriter(ostrm);
            }
            catch (Exception e)
            {
                Console.WriteLine("Cannot open RedirectV.txt for writing");
                Console.WriteLine(e.Message);
                return;
            }
            Console.SetOut(writer);

            // Read each scan in the RAW File
            for (int scan = firstScanNumber; scan <= lastScanNumber; scan++)
            {
                // Get the retention time for this scan number.  This is one of two comparable functions that will
                // convert between retention time and scan number.
                double time = rawFile.RetentionTimeFromScanNumber(scan);

                // Get the scan filter for this scan number
                var scanFilter = rawFile.GetFilterForScanNumber(scan);

                // Get the scan event for this scan number
                var scanEvent = rawFile.GetScanEventForScanNumber(scan);
               
                // Get the ionizationMode, MS2 precursor mass, collision energy, and isolation width for each scan
                if (scanFilter.MSOrder == MSOrderType.Ms2)
                {
                    // Get the reaction information for the first precursor
                    var reaction = scanEvent.GetReaction(0);

                    double precursorMass = reaction.PrecursorMass;
                    double collisionEnergy = reaction.CollisionEnergy;
                    double isolationWidth = reaction.IsolationWidth;
                    double monoisotopicMass = 0.0;
                    int masterScan = 0;
                    var ionizationMode = scanFilter.IonizationMode;
                    var order = scanFilter.MSOrder;

                    // Get the trailer extra data for this scan and then look for the monoisotopic m/z value in the 
                    // trailer extra data list
                    var trailerData = rawFile.GetTrailerExtraInformation(scan);

                    //Console.WriteLine(trailerData.Length);
                    //Console.WriteLine("is trailerData.Length and we have MS2 now");

                    for (int i = 0; i < trailerData.Length; i++)
                    {
                        if (trailerData.Labels[i] == "Monoisotopic M/Z:")
                        {
                            monoisotopicMass = Convert.ToDouble(trailerData.Values[i]);
                        }

                        if ((trailerData.Labels[i] == "Master Scan Number:") || (trailerData.Labels[i] == "Master Scan Number") || (trailerData.Labels[i] == "Master Index:"))
                        {
                            masterScan = Convert.ToInt32(trailerData.Values[i]);
                        }
                    }

                    if (outputData)
                    {
                        //Console.WriteLine(
                         //   "Scan number {0} @ time {1:F2} - Master scan = {2}, Ionization mode={3}, MS Order={4}, Precursor mass={5:F4}, Monoisotopic Mass = {6:F4}, Collision energy={7:F2}, Isolation width={8:F2} ",
                       //     scan, time, masterScan, ionizationMode, order, precursorMass, monoisotopicMass, collisionEnergy, isolationWidth);
                    }

                    var centroidStream = rawFile.GetCentroidStream(scan, false);

                    //get the time here
                    double timeInfo = rawFile.RetentionTimeFromScanNumber(scan);

                    //comment out the MS1 to do just the MS2
                    for (int i = 0; i < centroidStream.Intensities.Length; i++)
                    {
                    //we've added the MS1 outputs to the centroidStreamcentroid stream
                    //we'll also want to add the MS2 outputs to the centroid stream
                    //this is the ms2 stuff, now
                     //var stringPrint = centroidStream.Masses[i] + "\t" + centroidStream.Intensities[i] + "\t" + scan + "\t" + scanFilter.Polarity + "\t" + scanFilter.MSOrder + "\t" + reaction.PrecursorMass + "\t" + time + "\t" + timeInfo;

                     var stringPrint = centroidStream.Masses[i] + "\t" + centroidStream.Intensities[i]  + "\t" + scanFilter.Polarity + "\t" + scanFilter.MSOrder + "\t" + reaction.PrecursorMass + "\t"+ "\t" + timeInfo;

                     Console.WriteLine(stringPrint);
                    //Console.WriteLine(centroidStream.Masses[i]);
                    }


                }
                if (scanFilter.MSOrder == MSOrderType.Ms)
                {


                    

                    var centroidStream = rawFile.GetCentroidStream(scan, false);



                    double timeInfo = rawFile.RetentionTimeFromScanNumber(scan);
                    var scanInfo = rawFile.GetScanStatsForScanNumber(scan);


                    var scanProf = Scan.FromFile(rawFile, scan);



                    ////this is all the profile data we'll come back to, note well
                   // var segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scan, scanInfo);

                    //for (int i = 0; i < segmentedScan.Positions.Length; i++)
                   // {
                     //   Console.WriteLine("  {0} - {1:F4}, {2:F0}", i, segmentedScan.Positions[i], segmentedScan.Intensities[i]);
                    //}


                    //Console.WriteLine(scanInfo);
                    //Console.WriteLine("mass = {0}, intensity = {1}, centroid = ${2}, number = {3}", scanInfo.BasePeakMass,scanInfo.BasePeakIntensity, scanInfo.BasePeakIntensity, centroidStream.Masses.ToString(), scan);
                    //Console.WriteLine("[{0}]", string.Join(", ", centroidStream.Masses));
                    //Console.WriteLine("[{0}]", string.Join(", ", centroidStream.Intensities));


                    //Console.WriteLine( centroidStream.Intensities.Length);
                    //Console.WriteLine(centroidStream.Masses.Length);

                    //comment out the MS1 to do just the MS2 this is done on august 11th
                    for (int i = 0; i < centroidStream.Intensities.Length; i++)
                    {
                        //we've added the MS1 outputs to the centroidStreamcentroid stream
                        //we'll also want to add the MS2 outputs to the centroid stream

                        //timeInfo instea of scan - note well
                        var stringPrint = centroidStream.Masses[i] + "\t" + centroidStream.Intensities[i] + "\t" + timeInfo + "\t" + scanFilter.Polarity + "\t" + scanFilter.MSOrder;
                       Console.WriteLine(stringPrint);
                       // Console.WriteLine(centroidStream.Masses[i]);
                    }


    // Console.WriteLine(scanInfo.BasePeakMass);

   // var scanDependents = rawFile.GetScanDependents(scan, 5);

                    //Console.WriteLine(
                        //"Scan number {0} @ time {1:F2} - Instrument type={2}, Number dependent scans={3}",
                        //scan, time, scanDependents.RawFileInstrumentType, scanDependents.ScanDependentDetailArray.Length);
                        //Console.WriteLine(trailerData.Length);
                      //  Console.WriteLine("is trailerData.Length");
                }

                
                //var scanDependents = rawFile.GetScanDependents(scan, 5);

               // Console.WriteLine("Scan number {0} @ time {1:F2} - Instrument type={2}, Number dependent scans={3}", scan, time, scanDependents.RawFileInstrumentType, scanDependents.ScanDependentDetailArray.Length);
                


            }

            Console.SetOut(oldOut);
            writer.Close();
            ostrm.Close();
            Console.WriteLine("Done");

        }
    }
}
