/*
MIT License

Copyright (c) 2016 Thorsten Wagner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

package de.biomedical_imaging.ij.trajectory_classifier;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentFeature;

import de.biomedical_imaging.traJ.features.ActiveTransportParametersFeature;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;

public class ExportImportTools {
	
	public void exportTrajectoryDataAsCSV(ArrayList<? extends Trajectory> tracks, String path){
		String[] nextLine = null;
		try {
			CSVWriter writer = new CSVWriter(new FileWriter(path, false));
			nextLine = new String[]{"ID","X","Y","CLASS"};
			writer.writeNext(nextLine);
			
			for(int i = 0; i < tracks.size(); i++){
				Trajectory t = tracks.get(i);
				for(int j = 0; j < t.size(); j++){
					nextLine = new String[]{""+t.getID(),""+t.get(j).x,""+t.get(j).y,t.getType()};
					writer.writeNext(nextLine);
				}
			}
			writer.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void exportMSDDataAsCSV(ArrayList<? extends Trajectory> tracks, String path){
		String[] nextLine = null;
		try {
			double timelag = TraJClassifier_.getInstance().getTimelag();

			for(int j = 0; j < tracks.size(); j++){
				Trajectory t = tracks.get(j);

				int lastDotIndex = path.lastIndexOf('.');
				String s = path.substring(0, lastDotIndex ) + "-" + t.getType().replace('/', '_') + "-" + t.getID() + path.substring(lastDotIndex);

				CSVWriter writer = new CSVWriter(new FileWriter(s, false));

				int lagMin = 1;
				int lagMax = t.size()/3;

				double[] xData = new double[lagMax - lagMin + 1];
				double[] yData = new double[lagMax - lagMin + 1];
				double[] modelData = new double[lagMax - lagMin + 1];

				MeanSquaredDisplacmentFeature msdeval = new MeanSquaredDisplacmentFeature(t, lagMin);

				msdeval.setTrajectory(t);
				msdeval.setTimelag(lagMin);

				for(int i = lagMin; i < lagMax + 1; ++i) {
					msdeval.setTimelag(i);
					double msdhelp = msdeval.evaluate()[0];
					xData[i - lagMin] = i;
					yData[i - lagMin] = msdhelp;
				 }

				String modelString = "";

				if(t.getType().equals("SUBDIFFUSION")){
					PowerLawFeature pwf = new PowerLawFeature(t, 1/timelag,1, t.size()/3);
					double[] res = pwf.evaluate();
					double a = res[0];
					double D = res[1];

					modelString = "y=4*D*t^alpha";

					for(int i = lagMin; i < lagMax + 1; ++i) {
						modelData[i - lagMin] = 4.0 * D * Math.pow((double)i * timelag, a);
					 }
				} else if(t.getType().equals("CONFINED")){
					ConfinedDiffusionParametersFeature cfeature = new ConfinedDiffusionParametersFeature(t, timelag,TraJClassifier_.getInstance().isUseReducedModelConfinedMotion());
					double[] res= cfeature.evaluate();
					double a;
					double b;
					double c;
					double d;

					if(TraJClassifier_.getInstance().isUseReducedModelConfinedMotion()){
						a = res[0];
						b = 1;
						c = 1;
						d = res[1];

					} else{
						a = res[0];
						b = res[2];
						c = res[3];
						d = res[1];
					}

					if (Math.abs(1.0 - b) < 1.0E-5 && Math.abs(1.0 - a) < 1.0E-5) {
						modelString = "y=a*(1-exp(-4*D*t/a))";
					 } else {
						modelString = "y=a*(1-b*exp(-4*c*D*t/a))";
					 }

					for(int i = lagMin; i < lagMax + 1; ++i) {
						modelData[i - lagMin] = a * (1.0 - b * Math.exp(-4.0 * d * ((double)i * timelag / a) * c));
					}
				} else if(t.getType().equals("NORM. DIFFUSION")){
					RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1/timelag, 1, t.size()/3);
					double[] res =regest.evaluate();
					double diffusionCoefficient = res[0];
					double intercept = res[2];

					modelString = "y=4*D*t + a";

					for(int i = lagMin; i < lagMax + 1; ++i) {
						modelData[i - lagMin] = intercept + 4.0 * diffusionCoefficient * (double)i * timelag;
					}
				} else if(t.getType().equals("DIRECTED/ACTIVE")){
					ActiveTransportParametersFeature apf = new ActiveTransportParametersFeature(t, timelag);
					double[] ares = apf.evaluate();
					double diffusionCoefficient = ares[0];
					double velocity = ares[1];
					modelString = "y=4*D*t + (v*t)^2";

					for(int i = lagMin; i < lagMax + 1; ++i) {
						modelData[i - lagMin] = Math.pow(velocity * (i * timelag), 2) + 4 * diffusionCoefficient*(i * timelag);
					 }
				}

				if (modelString != "") {
					nextLine = new String[]{"LAG","MSD", "model " + modelString};
				} else {
					nextLine = new String[]{"LAG","MSD"};
				}
				writer.writeNext(nextLine);

				for (int k = lagMin; k < lagMax + 1; k++) {
					msdeval.setTimelag(k);
					double msdhelp = msdeval.evaluate()[0];
					if(modelString != "") {
						nextLine = new String[]{""+k,""+msdhelp,""+modelData[k - lagMin]};
					} else {
						nextLine = new String[]{""+k,""+msdhelp};
					}
					writer.writeNext(nextLine);
				}
				writer.close();
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public ArrayList<Trajectory> importTrajectoryDataFromCSV(String path){
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		try {
			CSVReader reader = new CSVReader(new FileReader(path));
			String[] nextLine;
			reader.readNext(); //READ HEADER!
			Trajectory t =null;
			int lastID = -1;
			while ((nextLine = reader.readNext()) != null) {
				int nextID = Integer.parseInt(nextLine[0]) ;
				double nextX = Double.parseDouble(nextLine[1]);
				double nextY = Double.parseDouble(nextLine[2]);
				String nextClass = nextLine[3];
				if(nextID==lastID){
					System.out.println();
					t.add(nextX, nextY, 0);
					lastID=nextID;
				}else{
					if(t!=null){
						tracks.add(t);
					}
					t = new Trajectory(2);
					t.setID(nextID);
					t.setType(nextClass);
					t.add(nextX, nextY, 0);
					lastID = nextID;
				}
		    }
			tracks.add(t);
			reader.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return tracks;
		
	}

}
