import casefilepy
import datfilepy
import CRN_Gen
#import ReactorGen_min_internalmfc_outletvalve_newton_odeint_energy
import PostProc
import os

from glob import glob
from shutil import copyfile


def read_case():
    casefilepy.convert("CFD_CRN.cas")


def read_data(symm_axis):
    datfilepy.readdata("CFD_CRN.dat", symm_axis)


def generate_crn(choice_react, criteria, derived, extra, addon, startpt, chem_mech, tol_lim, id_dict, data_num,
                 symm_axis, avg_dict):
    path = CRN_Gen.generate(choice_react, criteria, derived, extra, addon, startpt, chem_mech, tol_lim, id_dict,
                            data_num, symm_axis, avg_dict)
    return path


def solve_crn(path, choice_energy, chem_mech, id_dict, data_num, symm_axis, Peclet):
    ReactorGen_min_internalmfc_outletvalve_newton_odeint_energy.Gen(path, choice_energy, chem_mech, id_dict, data_num,
                                                                    symm_axis, Peclet)


def doe_reactor_num(criteria, chem_mech, tol_lim, id_dict, addon, num, startpt, data_num, symm_axis, choice_energy,
                    Peclet, avg_dict={}):
    for i in num:
        path = generate_crn(i, criteria, '', 0, addon, startpt, chem_mech, tol_lim, id_dict, data_num, symm_axis,
                            avg_dict)
        solve_crn(path, choice_energy, chem_mech, id_dict, data_num, symm_axis, Peclet)


def doe_startpt(num, criteria, addon, chem_mech, tol_lim, id_dict, startpt, data_num, symm_axis, variables,
                axial_locations, Peclet):
    for i in startpt:
        path = generate_crn(num, criteria, 'None', 0, addon, i, chem_mech, tol_lim, id_dict, data_num, symm_axis)
        solve_crn(path, 'off', chem_mech, id_dict, data_num, symm_axis, Peclet)
        diction = PostProc.Post(path, data_num, symm_axis)
        for variable in variables:
            PostProc.axial_plot(path, variable, diction, axial_locations, symm_axis, ['CRN'])


def write_results():
    loc = os.getcwd()
    try:
        os.mkdir(loc+'/Results')
    except:
        pass
    paths = glob('*/')
    print paths
    for p in paths:
        if p != 'Results\\':
            try:
                os.mkdir('Results/'+p)
            except:
                pass
            for files in glob(p+'*.vtu'):
                filename = files[len(p):]
                copyfile(files, 'Results/'+p+'/'+filename)
            for files in glob(p+'*.xlsx'):
                filename = files[len(p):]
                copyfile(files, 'Results/'+p+'/'+filename)
            for files in glob(p+'*readme.txt'):
                filename = files[len(p):]
                copyfile(files, 'Results/'+p+'/'+filename)


def main():
    # SET UP OF PLOTTING
    # Dictionary containing all the folders with data generated previously
    folder_root = 'D:/mddewit/Thesis/AGNES/SC_SFD_2D_'
    folder_root = os.getcwd()
    CRN_folder_dict = {
        folder_root + '153/1200_1955_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 1955 Monaghan',
        folder_root + '153/1200_1715_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 1715 Monaghan',
        folder_root + '153/1200_1559_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghann':
            'CRN PV1 $C\epsilon_1=1.53$ 1559 Monaghan',
        folder_root + '153/600_599_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 599 Monaghan',
        folder_root + '153/1050_1035_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 1035 Monaghan',
        folder_root + '153/250_249_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 249 Monaghan',
        folder_root + '153/125_151_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 151 Monaghan',
        folder_root + '153/65_141_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monaghan':
            'CRN PV1 $C\epsilon_1=1.53$ 141 Monaghan',
        folder_root + '153/141_141_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monen':
            'CRN PV1 $C\epsilon_1=1.53$ 141 Mon en',
        folder_root + '153/250_249_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monen':
            'CRN PV1 $C\epsilon_1=1.53$ 249 mon en',
        folder_root + '153/600_599_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monen':
            'CRN PV1 $C\epsilon_1=1.53$ 599 mon en',
        folder_root + '153/1036_1035_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monen':
            'CRN PV1 $C\epsilon_1=1.53$ 1035 mon en',
        folder_root + '153/1956_1955_Static Temperature_X-Coordinate_Mean Mixture Fraction__Monen':
            'CRN PV1 $C\epsilon_1=1.53$ 1955 mon en'}

    # Information about which folders should be used for plotting and where to save the plots
    plot_all = False
    folders_2_plot = []
    save_path = folder_root + '153/Figures/test'
    #path = folder_root + '/SandiaFlameD'
    #path = path.replace(os.sep, '/')

    # variables and axial locations to be plotted
    variables = ['Temperature[K]', 'MassfractionofCO2', 'MassfractionofCO', 'MassfractionofH2O', 'MassfractionofNO',
                 'MassfractionofO2', 'MassfractionofCH4']
    axial_locations = [0.0072, 0.0144, 0.0216, 0.108, 0.216, 0.324, 0.432, 0.54]

    # SET UP OF TEST CASE
    # chemical reaction mechanism to be used
    chem_mech = 'gri30.cti'

    # test case specific inputs from CFD
    startpt = 'fuel'
    data_num = 15
    symm_axis = 'X-Coordinate'
    id_dict = {'air': [21, {'O2': 0.23, 'N2': 0.77}, 291],
               'fuel': [17, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}, 294],
               'pilot': [20, {'O2': 0.059, 'N2': 0.735, 'H2O': 0.091, 'CO2': 0.111, 'OH': 0.001, 'NO': 0.003}, 1880],
               'top': [23],
               'side': [24]}

    """id_dict = {21: ['air', {'O2': 0.23, 'N2': 0.77}, 291],
               17: ['fuel', {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}, 294],
               20: ['pilot', {'O2': 0.059, 'N2': 0.735, 'H2O': 0.091, 'CO2': 0.111, 'OH': 0.001, 'NO': 0.003}, 1880],
               23: ['top'],
               24: ['side']}"""

    # tolerance and addition to the folder name to be used
    tol_lim = 100
    addon_list = ['']

    # START OF RUN
    # reading case and data information from CFD (comment out if already run for this specific CFD case and data)
    #read_case()
    #read_data(symm_axis)

    # information about the CRNs to be created (number of reactors, clustering criteria, use of energy equation)
    num_list = [1000]
    criteria_list = []
    for item in range(len(num_list)):
        criteria_list.append(['Static Temperature','Mass fraction of o2'])#, 'X-Coordinate', 'Mean Mixture Fraction'])
    energy_list = ['off']
    Peclet_list = [True]

    # dictionary containing zoning limits for the listed criteria (tuple) and the allowed difference in value for the
    # criteria in the reactor.
    # 'max'/'min' means the zone limit reaches the upper/lower bound of that variable,
    # 'none' means there is no clustering with respect to that criteria in that zone.
    # 'tol' means the tolerance specified above is used
    # 'full' is used to apply a specified tolerance. The tolerance is specified in the string with '*' between.
    # A string of a function can also be entered using 'mu' to apply scaling tolerances.
    # if avg_dict is not entered or left empty the clustering is done using only the tolerance supplied.
    for item in range(len(num_list)):
        avg_dict = {(0.0, 'max', 'min', 'max', 0.9, 1.0): ['none', 0.01, 0.01],
                    (0.0, 'max', 'min', 'max', 0.01, 0.1): [100.0, 'none', 'none'],
                    (0.0, 1800, 'min', 'max', 0.1, 0.9): [100.0, 'none', 'none'],
                    (1800, 2000, 'min', 'max', 0.1, 0.9): [50.0, 'none', 'none'],
                    (2000, 'max', 'min', 'max', 0.1, 0.9): [2.0, 'none', 'none'],
                    (0.0, 'max', 'min', 'max', 'min', 0.01): ['none', 0.2, 'none']}

        # avg_dict = {(0.0, 'max', 'min', 'max', 0.9, 1.0): ['none', '0.01*full', 0.01],
        #             (0.0, 'max', 'min', 'max', 0.01, 0.1): ['0.06*full', 'none', 'none'],
        #             (0.0, '0.9*max', 'min', 'max', 0.1, 0.9): ['0.06*full', 'none', 'none'],
        #             ('0.9*max', '0.95*max', 'min', 'max', 0.1, 0.9): ['0.03*full', 'none', 'none'],
        #             ('0.95*max', 'max', 'min', 'max', 0.1, 0.9): ['0.001*full', 'none', 'none'],
        #             (0.0, 'max', 'min', 'max', 'min', 0.01): ['none', '0.3*full', 'none']}
        #
        # avg_dict = {('min', 'max', 'min', 'max', 'min', 'max'): ['(-0.059*(mu[0]-cri_min[0])/avg_base[0]+0.06)*avg_base[0]',
        #                                                          '(-0.29*(mu[1]-cri_min[2])/avg_base[2]+0.3)*avg_base[1]',
        #                                                          '(-0.99*(mu[2]-cri_min[2])/avg_base[2]+1)*avg_base[2]']}

        # looping through run cases generating the CRN and solving it (comment out doe_reactor_num
        num = num_list[item]
        criteria = criteria_list[item]
        choice_energy = energy_list[item]
        addon = addon_list[item]
        Peclet = Peclet_list[item]
        avg_dict ={}
        doe_reactor_num(criteria, chem_mech, tol_lim, id_dict, addon, [num], startpt, data_num, symm_axis,
                         choice_energy, Peclet, avg_dict)
        #solve_crn(path, choice_energy, chem_mech, id_dict, data_num, symm_axis, Peclet)

    print "Writing Results"

    # START OF PLOTTING
    # create dictionaries of data from the test cases to be plotted and ParaView plots
    dictions = []
    result_legend = []
    result_legend2 = []

    if plot_all:
        for CRN_folder in CRN_folder_dict.keys():
            dictions.append(PostProc.Post(CRN_folder, data_num, symm_axis))
            result_legend.append(CRN_folder_dict[CRN_folder])
            result_legend2.append(CRN_folder_dict[CRN_folder])
            result_legend2.append(CRN_folder_dict[CRN_folder] + ' Volume average')

    else:
        for CRN_folder in folders_2_plot:
            dictions.append(PostProc.Post(CRN_folder, data_num, symm_axis))
            result_legend.append(CRN_folder_dict[CRN_folder])
            result_legend2.append(CRN_folder_dict[CRN_folder])
            result_legend2.append(CRN_folder_dict[CRN_folder] + ' Volume average')

    # generate case specific plots
    for variable in variables:
        try:
            vmin, vmax = PostProc.parallel_contour(save_path, variable, dictions, result_legend,
                                                   CRN_folder_dict[folders_2_plot[0]], axial_locations)
            PostProc.cross_contour(save_path, variable, dictions, axial_locations, symm_axis, result_legend, vmin, vmax,
                                   'CRN PV1 $C\epsilon_1=1.53$')
            integrals = PostProc.integral_plot('/'.join(save_path.split('/')[-2:]) + '/', variable, dictions,
                                               axial_locations, symm_axis, result_legend, 0.0072)
            PostProc.axial_plot(save_path, variable, dictions, axial_locations, symm_axis, result_legend2, integrals)
        except IndexError:
            continue

    print "Program Terminated Successfully"

    return 0


if __name__ == "__main__":
    main()
