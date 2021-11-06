% Generate all figures and tables for 
% "How the value of the environment controls persistence in visual search"
% by Traner, Bromberg-Martin, and Monosov (PLoS Computational Biology)

disp('Figure 1');
Fig1_upload;
drawnow;

disp('Figure 2');
Fig2_upload_general;
drawnow;

disp('Figure 3');
Fig3_upload_general;
drawnow;

disp('Figure 4');
Fig4_upload_general;
drawnow;

disp('Figure 5 and supp Fig 3');
Fig5_and_suppFigure3_general;
drawnow;

disp('supp Fig 1');
suppFigure1_general;
drawnow;

disp('supp Fig 2');
suppFigure2_general;
drawnow;

disp('supp Fig 4');
suppFigure4_general;
drawnow;

disp('supp Fig 5');
suppFigure5_general;
drawnow;

disp('Table');
tableGenerator;
