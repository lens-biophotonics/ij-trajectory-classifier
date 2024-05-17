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

import java.awt.Frame;
import java.awt.Menu;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.knowm.xchart.Chart;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.VisualizationUtils;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.ActiveTransportParametersFeature;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import ij.IJ;
import ij.WindowManager;
import ij.measure.ResultsTable;
import ij.text.TextPanel;

public class TraJResultsTable extends ResultsTable {
	
	boolean isParentTable;
	public TraJResultsTable() {
		this.isParentTable = false;
	}
	
	public TraJResultsTable(boolean isParentTable) {
		this.isParentTable = isParentTable;
	}
	
	@Override
	public void show(String windowTitle) {
		// TODO Auto-generated method stub
		super.show(windowTitle);
		final ResultsTable table = this;
		
		
		/*
		 * Plot trajectory
		 */
		MenuItem plotSelectedTrajectory = new MenuItem("Plot trajectory");
		plotSelectedTrajectory.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {

				Frame f = (Frame)WindowManager.getActiveWindow();
				if (f.getComponent(0) instanceof TextPanel){
					TextPanel p = (TextPanel) f.getComponent(0);
				
					int selectionStart = p.getSelectionStart();
					int selectionEnd = p.getSelectionEnd();
					ArrayList<Chart> charts = new ArrayList<Chart>();
					if(selectionStart>=0 && selectionStart==selectionEnd){
						int id = (int) table.getValue("ID", selectionStart);
						
						ArrayList<? extends Trajectory> cTracks;
						if(isParentTable){
							cTracks = TraJClassifier_.getInstance().getParentTrajectories();
						}else{
							cTracks = TraJClassifier_.getInstance().getClassifiedTrajectories();
						}
						Trajectory t = TrajectoryUtil.getTrajectoryByID(cTracks, id);
						Chart c = VisualizationUtils.getTrajectoryChart("Trajectory with ID " + id,t);
						charts.add(c);
						double timelag = TraJClassifier_.getInstance().getTimelag();
						if(t.getType().equals("SUBDIFFUSION")){
							PowerLawFeature pwf = new PowerLawFeature(t, 1/timelag,1, t.size()/3);
							double[] res = pwf.evaluate();
							c = VisualizationUtils.getMSDLineWithPowerModelChart(t, 1, t.size()/3, timelag, res[0], res[1]);
							charts.add(c);
							

						}else if(t.getType().equals("CONFINED")){
							ConfinedDiffusionParametersFeature cfeature = new ConfinedDiffusionParametersFeature(t, timelag,TraJClassifier_.getInstance().isUseReducedModelConfinedMotion());
							double[] res= cfeature.evaluate();
							if(TraJClassifier_.getInstance().isUseReducedModelConfinedMotion()){
								c = VisualizationUtils.getMSDLineWithConfinedModelChart(t, 1, t.size()/3, timelag, res[0], 1, 1, res[1]);
								charts.add(c);
							}else{
								c = VisualizationUtils.getMSDLineWithConfinedModelChart(t, 1, t.size()/3, timelag, res[0], res[2], res[3], res[1]);
								charts.add(c);
							}
						}else if(t.getType().equals("NORM. DIFFUSION")){
							RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1/timelag, 1, t.size()/3);
							double[] res =regest.evaluate();
							double dc = res[0];
							double intercept = res[2];
							c = VisualizationUtils.getMSDLineWithFreeModelChart(t, 1, t.size()/3, timelag, dc, intercept);
							charts.add(c);
						}else if(t.getType().equals("DIRECTED/ACTIVE")){
							ActiveTransportParametersFeature apf = new ActiveTransportParametersFeature(t, timelag);
							double[] ares = apf.evaluate();
							c = VisualizationUtils.getMSDLineWithActiveTransportModelChart(t, 1, t.size()/3, timelag, ares[0],ares[1]);
							charts.add(c);
						}else{
							c = VisualizationUtils.getMSDLineChart(t, 1, t.size()/3);
							charts.add(c);
						}
						VisualizationUtils.plotCharts(charts);
					}else if( selectionStart!=selectionEnd){
						IJ.showMessage("Plot of multiple trajectories is not possible");
					}
					else{
						IJ.showMessage("No trajectory selected");
					}
				}
				
			}
		});
		
		/*
		 * Export trajectory(s)
		 */
		MenuItem exportTrajectories = new MenuItem("Export trajetorie(s)");
		
		exportTrajectories.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				Frame f = (Frame)WindowManager.getActiveWindow();
				
				if (f.getComponent(0) instanceof TextPanel){
					TextPanel p = (TextPanel)f.getComponent(0);
				
					int selectionStart = p.getSelectionStart();
					int selectionEnd = p.getSelectionEnd();
					if(selectionStart == -1 && selectionEnd == -1){
						selectionStart = 0;
						selectionEnd = p.getResultsTable().getCounter();
					}
				
					ArrayList<Trajectory> selectedTrajectories = new ArrayList<Trajectory>();
					for( int i = selectionStart; i <= selectionEnd; i++){
						int id = (int) table.getValue("ID", i);
						
						ArrayList<? extends Trajectory> cTracks = null;
						if(isParentTable){
							cTracks = TraJClassifier_.getInstance().getParentTrajectories();
						}else{
							cTracks = TraJClassifier_.getInstance().getClassifiedTrajectories();
						}
						Trajectory t = TrajectoryUtil.getTrajectoryByID(cTracks, id);
						selectedTrajectories.add(t);
					}
					
					JFileChooser chooser=new JFileChooser();
					chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
					chooser.addChoosableFileFilter(new FileNameExtensionFilter("Comma seperated value (.csv)", "csv"));
					chooser.setAcceptAllFileFilterUsed(false);
					int c =chooser.showSaveDialog(null);
					if(c == JFileChooser.APPROVE_OPTION){
						String path=chooser.getSelectedFile().getAbsolutePath();
						if(!path.substring(path.length()-3, path.length()).equals("csv")){
							path += ".csv";
						}
				
						ExportImportTools eit = new ExportImportTools();
						eit.exportTrajectoryDataAsCSV(selectedTrajectories, path);
					}
					
				}
			}
		});

		/*
		 * Export MSD(s)
		 */
		MenuItem exportMSD = new MenuItem("Export MSD(s)");

		exportMSD.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				Frame f = (Frame)WindowManager.getActiveWindow();
				
				if (f.getComponent(0) instanceof TextPanel){
					TextPanel p = (TextPanel)f.getComponent(0);
				
					int selectionStart = p.getSelectionStart();
					int selectionEnd = p.getSelectionEnd();
					if(selectionStart == -1 && selectionEnd == -1){
						selectionStart = 0;
						selectionEnd = p.getResultsTable().getCounter();
					}

					ArrayList<Trajectory> selectedTrajectories = new ArrayList<Trajectory>();
					for( int i = selectionStart; i <= selectionEnd; i++){
						int id = (int) table.getValue("ID", i);

						ArrayList<? extends Trajectory> cTracks = null;
						if(isParentTable){
							cTracks = TraJClassifier_.getInstance().getParentTrajectories();
						}else{
							cTracks = TraJClassifier_.getInstance().getClassifiedTrajectories();
						}
						Trajectory t = TrajectoryUtil.getTrajectoryByID(cTracks, id);
						selectedTrajectories.add(t);
					}

					JFileChooser chooser=new JFileChooser();
					chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
					chooser.addChoosableFileFilter(new FileNameExtensionFilter("Comma seperated value (.csv)", "csv"));
					chooser.setAcceptAllFileFilterUsed(false);
					int c =chooser.showSaveDialog(null);
					if(c == JFileChooser.APPROVE_OPTION){
						String path=chooser.getSelectedFile().getAbsolutePath();
						if(!path.substring(path.length()-3, path.length()).equals("csv")){
							path += ".csv";
						}
				
						ExportImportTools eit = new ExportImportTools();
						eit.exportMSDDataAsCSV(selectedTrajectories, path);
					}
				}
			}
		});
		
		Menu traJ = new Menu("Trajectory classifier");
		traJ.add(plotSelectedTrajectory);
		traJ.add(exportTrajectories);
		traJ.add(exportMSD);
				//ResultsTable.getResultsWindow().getMenuBar().
		//Hinzufügen von Export Funktionen
		WindowManager.getFrame(windowTitle).getMenuBar().add(traJ);
	}

}
