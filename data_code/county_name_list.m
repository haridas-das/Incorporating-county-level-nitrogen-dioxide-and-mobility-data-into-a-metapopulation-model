function [state, m, name_origin_list, time_start, time_end, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list()

%% variable parameters 

state = 40;

m = 77;

name_origin_list = {'Adair County', 'Alfalfa County', 'Atoka County', 'Beaver County', 'Beckham County', 'Blaine County', 'Bryan County', 'Caddo County', 'Canadian County', 'Carter County', 'Cherokee County', 'Choctaw County', 'Cimarron County', 'Cleveland County', 'Coal County', 'Comanche County', 'Cotton County', 'Craig County', 'Creek County', 'Custer County', 'Delaware County', 'Dewey County', 'Ellis County', 'Garfield County', 'Garvin County', 'Grady County', 'Grant County', 'Greer County', 'Harmon County', 'Harper County', 'Haskell County', 'Hughes County', 'Jackson County', 'Jefferson County', 'Johnston County', 'Kay County', 'Kingfisher County', 'Kiowa County', 'Latimer County', 'Le Flore County', 'Lincoln County', 'Logan County', 'Love County', 'Major County', 'Marshall County', 'Mayes County', 'McClain County', 'McCurtain County', 'McIntosh County', 'Murray County', 'Muskogee County', 'Noble County', 'Nowata County', 'Okfuskee County', 'Oklahoma County', 'Okmulgee County', 'Osage County', 'Ottawa County', 'Pawnee County', 'Payne County', 'Pittsburg County', 'Pontotoc County', 'Pottawatomie County', 'Pushmataha County', 'Roger Mills County', 'Rogers County', 'Seminole County', 'Sequoyah County', 'Stephens County', 'Texas County', 'Tillman County', 'Tulsa County', 'Wagoner County', 'Washington County', 'Washita County', 'Woods County', 'Woodward County'};

time_start = 1;

time_end = 88;

top_five_pop_flow_nodes = [55, 72, 14, 9, 16];

%% constant parameters which we may change 

time_step_size = 0.10;                                                     % can compute N_time_steps

zeta_vec = [5, 6, 7, 8, 9, 10];                                            % put the diferent zeta values 
             
zeta = zeta_vec(1) ; % Dimensionless multiplicative factor that modifies the infectivity at the heterogeneous node.

Nv = 3;

gamma = 0.1429;                                                            % 1/7 wekly recovered

beta_1_vec = linspace(0,0.1428,1);                                         % simulate the infection rates

phi_0 = 0.0068;                                                             % mean flux value in the first 81 days i.e. 0.006795

end