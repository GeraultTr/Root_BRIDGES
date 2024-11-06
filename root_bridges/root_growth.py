from rhizodep.root_growth import RootGrowthModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

from openalea.mtg import *
from openalea.mtg.traversal import post_order
from numpy import pi, sqrt
import numpy as np
import pandas as pd
import os


family = "growth"


@dataclass
class RootGrowthModelCoupled(RootGrowthModel):
    """

    NOTE : double names in methods are forbiden as they will be overwritten in Choregrapher's resolution
    However, it is a behavios of interest here when inheriting the class to edit it.
    """
    # INPUTS
    
    nitrate_transporters_affinity_factor: float = declare(default=1., unit="mol.s-1", unit_comment="of nitrates", description="nitrate_transporters_affinity_factor, introduced to account for NRT1 signalling function when going through LATS regime", 
                                                    min_value="", max_value="", value_comment="", references="Remans et al 2006", DOI="", 
                                                    variable_type="input", by="model_nitrogen", state_variable_type="intensive", edit_by="user")

    # STATE VARIABLES
    amino_acids_consumption_by_growth: float = declare(default=0., unit="mol.s-1", unit_comment="", description="amino_acids consumption rate by growth processes", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                                    variable_type="state_variable", by="model_growth", state_variable_type="extensive", edit_by="user")
    
    # PARAMETERS
    Km_elongation_amino_acids: float = declare(default=1250 * 1e-6 / 5, unit="mol.g-1", unit_comment="of amino_acids", description="Affinity constant for root elongation regarding amino_acids consumption",
                                                    min_value="", max_value="", value_comment="TODO : actualize", references="According to Barillot et al. (2016b): Km for root growth is 1250 umol C g-1 for sucrose. According to Gauthier et al (2020): Km for regulation of the RER by sucrose concentration in hz = 100-150 umol C g-1", DOI="",
                                                    variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")
    Km_nodule_thickening_amino_acids: float = declare(default=1250 * 1e-6 / 5. * 100, unit="mol.g-1", unit_comment="of amino_acids", description="Affinity constant for nodule thickening regarding amino_acids consumption", 
                                                    min_value="", max_value="", value_comment="Km_elongation * 100, TODO : actualize", references="", DOI="",
                                                    variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")
    struct_mass_N_content: float = declare(default=0.03 / 14, unit="mol.g-1", unit_comment="of organic nitrogen", description="organic nitrogen content of structural mass",
                                                    min_value="", max_value="", value_comment="Increased from 3% to consider active protein requirements in structural mass", references="We assume that the structural mass contains 3% of N. (??)", DOI="",
                                                    variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")
    yield_growth_N: float = declare(default=1., unit="adim", unit_comment="mol of N per mol of C used for structural mass", description="Growth yield of amino acids", 
                                                    min_value="", max_value="", value_comment="No deamination considered during biosynthesis, so we stick to 1.", references="", DOI="",
                                                    variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")
    r_Nm_AA: float =     declare(default=1.4, unit="adim", unit_comment="mol of nitrogen per mol of amino acids", description="concentration stoechiometric ratio between mineral nitrogen and amino acids in roots", 
                                min_value="", max_value="", value_comment="TODO : check estimation", references="", DOI="",
                                variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")
    r_C_AA: float =     declare(default=5, unit="adim", unit_comment="mol of carbon per mol of amino acids", description="concentration stoechiometric ratio between carbon and amino acids in roots", 
                                min_value="", max_value="", value_comment="Based on glutamic acid", references="", DOI="",
                                variable_type="parameter", by="model_growth", state_variable_type="", edit_by="user")


    def __init__(self, g=None ,time_step=3600, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)


    def post_growth_updating(self):
        self.vertices = self.g.vertices(scale=self.g.max_scale())
        for vid in self.vertices:
            if vid not in self.amino_acids_consumption_by_growth.keys():
                parent = self.g.parent(vid)
                # we partition the initial flow in the parent accounting for mass fraction
                # We use struct_mass, the resulting structural mass after growth
                mass_fraction = self.struct_mass[vid] / (self.struct_mass[vid] + self.struct_mass[parent])
                self.amino_acids_consumption_by_growth.update({vid: self.amino_acids_consumption_by_growth[parent] * mass_fraction,
                                            parent: self.amino_acids_consumption_by_growth[parent] * (1 - mass_fraction)})

    # SUBDIVISIONS OF THE SCHEDULING LOOP
    # -----------------------------------
    @stepinit
    def reinitializing_growth_variables(self):
        """
        This function re-initializes different growth-related variables (e.g. potential growth variables).
        EDIT : Added amino acids growth and elongation variables for reinitialization TODO : Really usefull for all? in Some funcs, it seems repeated.
        :return:
        """
        # We cover all the vertices in the MTG:
        for vid in self.g.vertices_iter(scale=1):
            # n represents the vertex:
            n = self.g.node(vid)

            # We set to 0 the growth-related variables:
            n.hexose_consumption_by_growth_amount = 0.
            n.hexose_consumption_by_growth = 0.
            n.hexose_possibly_required_for_elongation = 0.
            n.amino_acids_consumption_by_growth_amount = 0.
            n.amino_acids_consumption_by_growth = 0.
            n.amino_acids_possibly_required_for_elongation = 0.
            n.resp_growth = 0.
            n.struct_mass_produced = 0.
            n.root_hairs_struct_mass_produced = 0.
            n.hexose_growth_demand = 0.
            n.amino_acids_growth_demand = 0.
            n.actual_elongation = 0.
            n.actual_elongation_rate = 0.

            n.hexose_consumption_by_fungus = 0.

            # We make sure that the initial values of length, radius and struct_mass are correctly initialized:
            n.initial_length = n.length
            n.initial_radius = n.radius
            n.potential_radius = n.radius
            n.theoretical_radius = n.radius
            n.initial_struct_mass = n.struct_mass
            n.initial_living_root_hairs_struct_mass = n.living_root_hairs_struct_mass
        return

    # Function for calculating root elongation:
    def elongated_length(self, element, initial_length: float, radius: float, C_hexose_root: float, elongation_time_in_seconds: float):
        """
        This function computes a new length (m) based on the elongation process described by ArchiSimple and regulated by
        the available concentration of hexose.
        It has been modified in Root-BRIDGES to introduce a regulation by 
        :param initial_length: the initial length (m)
        :param radius: radius (m)
        :param C_hexose_root: the concentration of hexose available for elongation (mol of hexose per gram of strctural mass)
        :param elongation_time_in_seconds: the period of elongation (s)
        :return: the new elongated length
        """

        # If we keep the classical ArchiSimple rule:
        if self.simple_growth_duration:
            # Then the elongation is calculated following the rules of Pages et al. (2014):
            elongation = self.EL * 2. * radius * elongation_time_in_seconds
        else:
            # Otherwise, we additionally consider a limitation of the elongation according to the local concentration of hexose,
            # based on a Michaelis-Menten formalism:
            if C_hexose_root > self.C_hexose_min_for_elongation and element.AA > 0:
                # michaelis_menten_limitation = ((1 + self.Km_elongation) / C_hexose_root) * ((1 + self.Km_elongation_amino_acids) / element.AA)
                michaelis_menten_limitation = ((C_hexose_root / (C_hexose_root + self.Km_elongation)) + (element.AA / (element.AA + self.Km_elongation_amino_acids))) / 2
                #print("MM", michaelis_menten_limitation)
                potential_elongation = self.EL * 2. * radius * elongation_time_in_seconds
                elongation = potential_elongation * michaelis_menten_limitation
            else:
                print(f"For element {element.index()}, no elongation, negative concentrations!! ", C_hexose_root, element.AA)
                elongation = 0.
        
        # We calculate the new potential length corresponding to this elongation:
        new_length = initial_length + elongation
        if new_length < initial_length:
            print("!!! ERROR: There is a problem of elongation, with the initial length", initial_length,
                  " and the radius", radius, "and the elongation time", elongation_time_in_seconds)
        return new_length

    # Function for calculating the amount of C to be used in neighbouring elements for sustaining root elongation:
    def calculating_supply_for_elongation(self, element):
        """
        This function computes the list of root elements that can supply C as hexose for sustaining the elongation
        of a given element, as well as their structural mass and their amount of available hexose.
        EDIT : Sustaining the need for N when a root apex should elongate has been added.
        :param element: the element for which we calculate the possible supply of C for its elongation
        :return: three lists containing the indices of elements, their hexose amount (mol of hexose) and their structural mass (g).
        """

        n = element

        # We initialize each amount of hexose available for growth:
        n.hexose_possibly_required_for_elongation = 0.
        n.amino_acids_possibly_required_for_elongation = 0.
        n.struct_mass_contributing_to_elongation = 0.

        # We initialize empty lists:
        list_of_elongation_supporting_elements = []
        list_of_elongation_supporting_elements_hexose = []
        list_of_elongation_supporting_elements_amino_acids = []
        list_of_elongation_supporting_elements_mass = []

        # We then calculate the length of an apical zone of a fixed length which can provide the amount of hexose required for growth:
        growing_zone_length = self.growing_zone_factor * n.radius
        # We calculate the corresponding volume to which this length should correspond based on the diameter of this apex:
        supplying_volume = growing_zone_length * n.radius ** 2 * pi

        # We start counting the hexose at the apex:
        index = n.index()
        current_element = n

        # We initialize a temporary variable that will be used as a counter:
        remaining_volume = supplying_volume

        # As long the remaining volume is not zero:
        while remaining_volume > 0:

            # If the volume of the current element is lower than the remaining volume:
            if remaining_volume > current_element.volume:
                # We make sure to include in the list of supplying elements only elements with a positive length
                # (e.g. NOT the elements of length 0 that support seminal or adventitious roots):
                if current_element.length > 0.:
                    # We add to the amount of hexose available all the hexose in the current element
                    # (EXCLUDING sugars in the living root hairs):
                    # TODO: Should the C from root hairs be used for helping roots to grow?
                    hexose_contribution = current_element.C_hexose_root * current_element.struct_mass
                    amino_acids_contribution = current_element.AA * current_element.struct_mass
                    n.hexose_possibly_required_for_elongation += hexose_contribution
                    n.amino_acids_possibly_required_for_elongation += amino_acids_contribution
                    n.struct_mass_contributing_to_elongation += current_element.struct_mass
                    # We record the index of the contributing element:
                    list_of_elongation_supporting_elements.append(index)
                    # We record the amount of hexose that the current element can provide:
                    list_of_elongation_supporting_elements_hexose.append(hexose_contribution)
                    # We record the amount of amino acids that the current element can provide:
                    list_of_elongation_supporting_elements_amino_acids.append(amino_acids_contribution)
                    # We record the structural mass from which the current element contributes:
                    list_of_elongation_supporting_elements_mass.append(current_element.struct_mass)
                    # We subtract the volume of the current element to the remaining volume:
                    remaining_volume = remaining_volume - current_element.volume

                # And we try to move the index to the segment preceding the current element:
                index_attempt = self.g.Father(index, EdgeType='<')
                # If there is no father element on this axis:
                if index_attempt is None:
                    # Then we try to move to the mother root, if any:
                    index_attempt = self.g.Father(index, EdgeType='+')
                    # If there is no such root:
                    if index_attempt is None:
                        # Then we exit the loop here:
                        break
                # We set the new index:
                index = index_attempt
                # We define the new element to consider according to the new index:
                current_element = self.g.node(index)
                
            # Otherwise, this is the last preceding element to consider:
            else:
                # We finally add to the amount of hexose available for elongation a part of the hexose of the current element:
                hexose_contribution = current_element.C_hexose_root * current_element.struct_mass \
                                      * remaining_volume / current_element.volume
                amino_acids_contribution = current_element.AA * current_element.struct_mass \
                                      * remaining_volume / current_element.volume
                n.hexose_possibly_required_for_elongation += hexose_contribution
                n.amino_acids_possibly_required_for_elongation += amino_acids_contribution
                n.struct_mass_contributing_to_elongation += current_element.struct_mass \
                                                            * remaining_volume / current_element.volume
                # We record the index of the contributing element:
                list_of_elongation_supporting_elements.append(index)
                # We record the amount of hexose that the current element can provide:
                list_of_elongation_supporting_elements_hexose.append(hexose_contribution)
                # We record the amount of amino acids that the current element can provide:
                list_of_elongation_supporting_elements_amino_acids.append(amino_acids_contribution)
                # We record the structural mass from which the current element contributes:
                list_of_elongation_supporting_elements_mass.append(
                    current_element.struct_mass * remaining_volume / current_element.volume)
                # And the remaining volume to consider is set to 0:
                remaining_volume = 0.
                # And we exit the loop here:
                break

        # We record the average concentration in hexose of the whole zone of hexose supply contributing to elongation:
        if n.struct_mass_contributing_to_elongation > 0.:
            n.growing_zone_C_hexose_root = n.hexose_possibly_required_for_elongation / n.struct_mass_contributing_to_elongation
        else:
            print("!!! ERROR: the mass contributing to elongation in element", n.index(), "of type", n.type, "is",
                  n.struct_mass_contributing_to_elongation,
                  "g, and its structural mass is", n.struct_mass, "g!")
            n.growing_zone_C_hexose_root = 0.

        n.list_of_elongation_supporting_elements = list_of_elongation_supporting_elements
        n.list_of_elongation_supporting_elements_hexose = list_of_elongation_supporting_elements_hexose
        n.list_of_elongation_supporting_elements_amino_acids = list_of_elongation_supporting_elements_amino_acids
        n.list_of_elongation_supporting_elements_mass = list_of_elongation_supporting_elements_mass


    def primordium_formation(self, apex, elongation_rate=0.):
        """
        This function considers the formation of a primordium on a root apex, and, if possible, creates this new element
        of length 0.
        :param apex: the apex on which the primordium may be created
        :param elongation_rate: the rate of elongation of the apex over the time step during which the primordium may be created (m s-1)
        :return: the new primordium
        """

        # NOTE: This function has to be called AFTER the actual elongation of the apex has been done and the distance
        # between the tip of the apex and the last ramification (dist_to_ramif) has been increased!

        # CALCULATING AN EQUIVALENT OF THERMAL TIME:
        # -------------------------------------------

        # We calculate a coefficient that will modify the different "ages" experienced by roots according to soil
        # temperature assuming a linear relationship (this is equivalent as the calculation of "growth degree-days):
        temperature_time_adjustment = self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=apex.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)

        # OPERATING PRIMORDIUM FORMATION:
        # --------------------------------
        # We initialize the new_apex that will be returned by the function:
        new_apex = []

        # VERIFICATION: We make sure that no lateral root has already been form on the present apex.
        # We calculate the number of children of the apex (it should be 0!):
        n_children = len(apex.children())
        # If there is at least one children, it means that there is already one lateral primordium or root on the apex:
        if n_children >= 1:
            # Then we don't add any primordium and simply return the unaltered apex:
            new_apex.append(apex)
            return new_apex

        # If the order of the current apex is too high, we forbid the formation of a new primordium of higher order:
        if self.root_order_limitation and apex.root_order + 1 > self.root_order_treshold:
            # Then we don't add any primordium and simply return the unaltered apex:
            new_apex.append(apex)
            return new_apex

        # We first calculate the radius that the primordium may have. This radius is drawn from a normal distribution
        # whose mean is the value of the mother root diameter multiplied by RMD, and whose standard deviation is
        # the product of this mean and the coefficient of variation CVDD (Pages et al. 2014).
        # We also set the root angles depending on random:
        if self.random:
            # The seed used to generate random values is defined according to a parameter random_choice and the index of the apex:
            np.random.seed(self.random_choice * apex.index())
            potential_radius = np.random.normal((apex.radius - self.Dmin / 2.) * self.RMD + self.Dmin / 2.,
                                                ((apex.radius - self.Dmin / 2.) * self.RMD + self.Dmin / 2.) * self.CVDD)
            apex_angle_roll = abs(np.random.normal(120, 10))
            if apex.root_order == 1:
                primordium_angle_down = abs(np.random.normal(45, 10))
            else:
                primordium_angle_down = abs(np.random.normal(70, 10))
            primordium_angle_roll = abs(np.random.normal(5, 5))
        else:
            potential_radius = (apex.radius - self.Dmin / 2) * self.RMD + self.Dmin / 2.
            apex_angle_roll = 120
            if apex.root_order == 1:
                primordium_angle_down = 45
            else:
                primordium_angle_down = 70
            primordium_angle_roll = 5

        # We add a condition, i.e. potential_radius should not be larger than the one from the mother root:
        if potential_radius > apex.radius:
            potential_radius = apex.radius

        # If the distance between the apex and the last emerged root is higher than the inter-primordia distance
        # AND if the potential radius is higher than the minimum diameter:
        if apex.dist_to_ramif > self.IPD and potential_radius >= self.Dmin and potential_radius <= apex.radius:
            # The distance that the tip of the apex has covered since the actual primordium formation is calculated:
            elongation_since_last_ramif = apex.dist_to_ramif - self.IPD
            # TODO: Is the third condition really relevant?

            # A specific rolling angle is attributed to the parent apex:
            apex.angle_roll = apex_angle_roll

            # We verify that the apex has actually elongated:
            if apex.actual_elongation > 0 and elongation_rate > 0.:
                # Then the actual time since the primordium must have been formed is precisely calculated
                # according to the actual growth of the parent apex since primordium formation,
                # taking into account the actual growth rate of the parent defined as
                # apex.actual_elongation / time_step_in_seconds
                actual_time_since_formation = elongation_since_last_ramif / elongation_rate
            else:
                actual_time_since_formation = 0.

            # And we add the primordium of a possible new lateral root:
            ramif = self.ADDING_A_CHILD(mother_element=apex, edge_type='+', label='Apex',
                                        type='Normal_root_before_emergence',
                                        root_order=apex.root_order + 1,
                                        angle_down=primordium_angle_down,
                                        angle_roll=primordium_angle_roll,
                                        length=0.,
                                        radius=potential_radius,
                                        identical_properties=False,
                                        nil_properties=True)
            if apex.nitrate_transporters_affinity_factor:
                lateral_elongation_possibility = max(apex.nitrate_transporters_affinity_factor * self.main_roots_growth_extender, 1)
            else:
                lateral_elongation_possibility = 1
            # We specifically recomputes the growth duration:
            ramif.growth_duration = self.GDs * (2. * ramif.radius) ** 2 * lateral_elongation_possibility
            # ramif.growth_duration = self.calculate_growth_duration(radius=ramif.radius, index=ramif.index(),
            #                                                        root_order=ramif.root_order)
            # We specify the exact time since formation:
            ramif.actual_time_since_primordium_formation = actual_time_since_formation
            ramif.thermal_time_since_primordium_formation = actual_time_since_formation * temperature_time_adjustment
            # And the new distance between the parent apex and the last ramification is redefined,
            # by taking into account the actual elongation of apex since the child formation:
            apex.dist_to_ramif = elongation_since_last_ramif
            # # We also put in memory the index of the child:
            # apex.lateral_primordium_index = ramif.index()
            # We add the apex and its ramif in the list of apices returned by the function:
            new_apex.append(apex)
            new_apex.append(ramif)

        return new_apex


    # Function calculating the potential development of a root segment:
    def potential_segment_development(self, segment):
        """
        This function considers a root segment, i.e. a root element that can thicken but not elongate, and calculates its
        potential increase in radius according to the pipe model (possibly regulated by C availability), and its possible death.

        EDIT : Added a regulation of root thickening with the availability of amino acids in the root segment through a bi-michaelian.

        :param segment: the segment to be considered
        :return: the updated segment
        """

        # We initialize an empty list that will contain the new segment to be returned:
        new_segment = []
        # We record the current radius and length prior to growth as the initial radius and length:
        segment.initial_radius = segment.radius
        segment.initial_length = segment.length
        # We initialize the properties "potential_radius" and "potential_length":
        segment.theoretical_radius = segment.radius
        segment.potential_radius = segment.radius
        segment.potential_length = segment.length

        # CASE 1: THE SEGMENT IS A NODULE:
        # ################################
        # NOTE: a nodule is considered here as a tumor which grows radially by feeding from root hexose, but does not produce
        # new root axes.

        if segment.type == "Root_nodule":
            # We consider the amount of hexose available in the nodule AND in the parent segment
            # (EXCLUDING the amount of hexose in living rot hairs):
            # TODO: Should the C from root hairs be used for helping nodules to grow?
            index_parent = self.g.Father(segment.index(), EdgeType='+')
            parent = self.g.node(index_parent)
            segment.hexose_available_for_thickening = parent.C_hexose_root * parent.struct_mass \
                                                      + segment.C_hexose_root * segment.struct_mass
            segment.amino_acids_available_for_thickening = parent.AA * parent.struct_mass \
                                                      + segment.AA * segment.struct_mass
            # We calculate an average concentration of hexose that will help to regulate nodule growth:
            C_hexose_regulating_nodule_growth = segment.hexose_available_for_thickening / (
                    parent.struct_mass + segment.struct_mass)
            N_amino_acids_regulating_nodule_growth = segment.amino_acids_available_for_thickening / (
                    parent.struct_mass + segment.struct_mass)
            # We modulate the relative increase in radius by the amount of C available in the nodule:
            thickening_rate = self.relative_nodule_thickening_rate_max / (
                ((1+self.Km_nodule_thickening) / C_hexose_regulating_nodule_growth)*((1+self.Km_nodule_thickening_amino_acids) / N_amino_acids_regulating_nodule_growth)
            )
            
            # We calculate a coefficient that will modify the rate of thickening according to soil temperature
            # assuming a linear relationship (this is equivalent as the calculation of "growth degree-days):
            thickening_rate = thickening_rate * self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=segment.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)
            segment.theoretical_radius = segment.radius * (1 + thickening_rate * self.time_step_in_seconds)
            if segment.theoretical_radius > self.nodule_max_radius:
                segment.potential_radius = self.nodule_max_radius
            else:
                segment.potential_radius = segment.theoretical_radius
            # We add the modified segment to the list of new segments, and we quit the function here:
            new_segment.append(segment)
            return new_segment

        # CASE 2: THE SEGMENT IS NOT A NODULE:
        ######################################

        # We initialize internal variables:
        son_section = 0.
        sum_of_lateral_sections = 0.
        number_of_actual_children = 0.
        death_count = 0.
        list_of_times_since_death = []

        # We define the amount of hexose available for thickening
        # (EXCLUDING the amount of hexose in living root hairs):
        # TODO: Should the C from root hairs be used for helping nodules to grow?
        segment.hexose_available_for_thickening = segment.C_hexose_root * segment.struct_mass
        segment.amino_acids_available_for_thickening = segment.AA * segment.struct_mass

        # CALCULATING AN EQUIVALENT OF THERMAL TIME:
        # ------------------------------------------

        # We calculate a coefficient that will modify the different "ages" experienced by roots according to soil
        # temperature assuming a linear relationship (this is equivalent as the calculation of "growth degree-days):
        temperature_time_adjustment = self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=segment.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)

        # CHECKING WHETHER THE APEX OF THE ROOT AXIS HAS STOPPED GROWING:
        # ---------------------------------------------------------------

        # We look at the apex of the axis to which the segment belongs (i.e. we get the last element of the axis):
        index_apex = self.g.Axis(segment.index())[-1]
        apex = self.g.node(index_apex)
        # print("For segment", segment.index(), "the terminal index is", index_apex, "and has the type", apex.label)
        # Depending on the type of the apex, we adjust the type of the segment on the same axis:
        if apex.type == "Just_stopped":
            segment.type = "Just_stopped"
        elif apex.type == "Stopped":
            segment.type = "Stopped"

        # CHECKING POSSIBLE ROOT SEGMENT DEATH:
        # -------------------------------------

        # For each child of the segment:
        for child in segment.children():

            # Then we add one child to the actual number of children:
            number_of_actual_children += 1

            if child.radius < 0. or child.potential_radius < 0.:
                print("!!! ERROR: the radius of the element", child.index(), "is negative!")
            # If the child belongs to the same axis:
            if child.edge_type == '<':
                # Then we record the THEORETICAL section of this child:
                son_section = child.theoretical_radius ** 2 * pi
                # # Then we record the section of this child:
                # son_section = child.radius * child.radius * pi
            # Otherwise if the child is the element of a lateral root AND if this lateral root has already emerged
            # AND the lateral element is not a nodule:
            elif child.edge_type == '+' and child.length > 0. and child.type != "Root_nodule":
                # We add the POTENTIAL section of this child to a sum of lateral sections:
                sum_of_lateral_sections += child.theoretical_radius ** 2 * pi
                # # We add the section of this child to a sum of lateral sections:
                # sum_of_lateral_sections += child.radius ** 2 * pi

            # If this child has just died or was already dead:
            if child.type == "Just_dead" or child.type == "Dead":
                # Then we add one dead child to the death count:
                death_count += 1
                # And we record the exact time since death:
                list_of_times_since_death.append(child.actual_time_since_death)

        # If each child in the list of children has been recognized as dead or just dead:
        if death_count == number_of_actual_children:
            # If the investigated segment was already declared dead at the previous time step:
            if segment.type == "Just_dead" or segment.type == "Dead":
                # Then we transform its status into "Dead"
                segment.type = "Dead"
            else:
                # Then the segment has to die:
                segment.type = "Just_dead"
        # Otherwise, at least one of the children axis is not dead, so the father segment should not be dead

        # REGULATION OF RADIAL GROWTH BY AVAILABLE CARBON:
        # ------------------------------------------------
        # If the radial growth is possible:
        if self.radial_growth:
            # The radius of the root segment is defined according to the pipe model.
            # In ArchiSimp9, the radius is increased by considering the sum of the sections of all the children,
            # by adding a fraction (SGC) of this sum of sections to the current section of the parent segment,
            # and by calculating the new radius that corresponds to this new section of the parent:
            segment.theoretical_radius = sqrt(son_section / pi + self.SGC * sum_of_lateral_sections / pi)
            # However, if the net difference is below 0.1% of the initial radius:
            if (segment.theoretical_radius - segment.initial_radius) <= 0.001 * segment.initial_radius:
                # Then the potential radius is set to the initial radius:
                segment.theoretical_radius = segment.initial_radius
            # If we consider simple ArchiSimple rules:
            if self.simple_growth_duration:
                # Then the potential radius to form is equal to the theoretical one determined by geometry:
                segment.potential_radius = segment.theoretical_radius
            # Otherwise, if we don't strictly follow simple ArchiSimple rules and if there can be an increase in radius:
            elif segment.length > 0. and segment.theoretical_radius > segment.radius and segment.C_hexose_root > self.C_hexose_min_for_thickening:
                # We calculate the maximal increase in radius that can be achieved over this time step,
                # based on a Michaelis-Menten formalism that regulates the maximal rate of increase
                # according to the amount of hexose available:
                thickening_rate = self.relative_root_thickening_rate_max / (
                    ((1+self.Km_nodule_thickening) / segment.C_hexose_root) + ((1+self.Km_nodule_thickening_amino_acids) / segment.AA)) / 2

                # We calculate a coefficient that will modify the rate of thickening according to soil temperature
                # assuming a linear relationship (this is equivalent as the calculation of "growth degree-days):
                thickening_rate = thickening_rate * self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=segment.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)
                # The maximal possible new radius according to this regulation is therefore:
                new_radius_max = (1 + thickening_rate * self.time_step_in_seconds) * segment.initial_radius
                # If the potential new radius is higher than the maximal new radius:
                if segment.theoretical_radius > new_radius_max:
                    # Then potential thickening is limited up to the maximal new radius:
                    segment.potential_radius = new_radius_max
                # Otherwise, the potential radius to achieve is equal to the theoretical one:
                else:
                    segment.potential_radius = segment.theoretical_radius
            # And if the segment corresponds to one of the elements of length 0 supporting one seminal or adventitious root:
            if segment.type == "Support_for_seminal_root" or segment.type == "Support_for_adventitious_root":
                # Then the radius is directly increased, as this element will not be considered in the function calculating actual growth:
                segment.radius = segment.potential_radius

        # UPDATING THE DIFFERENT TIMES:
        # ------------------------------

        # We increase the various time variables:
        segment.actual_time_since_primordium_formation += self.time_step_in_seconds
        segment.actual_time_since_emergence += self.time_step_in_seconds
        segment.actual_time_since_cells_formation += self.time_step_in_seconds
        segment.thermal_time_since_primordium_formation += self.time_step_in_seconds * temperature_time_adjustment
        segment.thermal_time_since_emergence += self.time_step_in_seconds * temperature_time_adjustment
        segment.thermal_time_since_cells_formation += self.time_step_in_seconds * temperature_time_adjustment

        if segment.type == "Just_stopped":
            segment.actual_time_since_growth_stopped = apex.actual_time_since_growth_stopped
            segment.thermal_time_since_growth_stopped = apex.actual_time_since_growth_stopped * temperature_time_adjustment
        if segment.type == "Stopped":
            segment.actual_time_since_growth_stopped += self.time_step_in_seconds
            segment.thermal_time_since_growth_stopped += self.time_step_in_seconds * temperature_time_adjustment
        if segment.type == "Just_dead":
            segment.actual_time_since_growth_stopped += self.time_step_in_seconds
            segment.thermal_time_since_growth_stopped += self.time_step_in_seconds * temperature_time_adjustment
            # AVOIDING PROBLEMS - We check that the list of times_since_death is not empty:
            if list_of_times_since_death:
                segment.actual_time_since_death = min(list_of_times_since_death)
            else:
                segment.actual_time_since_death = 0.
            segment.thermal_time_since_death = segment.actual_time_since_death * temperature_time_adjustment
        if segment.type == "Dead":
            segment.actual_time_since_growth_stopped += self.time_step_in_seconds
            segment.thermal_time_since_growth_stopped += self.time_step_in_seconds * temperature_time_adjustment
            segment.actual_time_since_death += self.time_step_in_seconds
            segment.thermal_time_since_death += self.time_step_in_seconds * temperature_time_adjustment

        new_segment.append(segment)
        return new_segment
    
    
    # Actual elongation, radial growth and growth respiration of root elements:
    @actual
    @state
    def actual_growth_and_corresponding_respiration(self):
        """
        This function defines how a segment, an apex and possibly an emerging root primordium will grow according to the amount
        of hexose present in the segment, taking into account growth respiration based on the model of Thornley and Cannell
        (2000). The calculation is based on the values of potential_radius, potential_length, lateral_root_emergence_possibility
        and emergence_cost defined in each element by the module "POTENTIAL GROWTH".
        The function returns the MTG "g" with modified values of radius and length of each element, the possibility of the
        emergence of lateral roots, and the cost of growth in terms of hexose consumption.

        EDIT : explicited a cost of amino acids associated to the actual growth, 
        and considered how to limit the actual growth when not enough amino acids are available through the limiting factor strategy between C and N.

        :return:
        """

        # TODO FOR TRISTAN: In a second step, you may consider how the emergence of primordia may depend on N availability in the root or in the soil.

        # PROCEEDING TO ACTUAL GROWTH:
        # -----------------------------

        # We have to cover each vertex from the apices up to the base one time:
        root_gen = self.g.component_roots_at_scale_iter(self.g.root, scale=1)
        root = next(root_gen)
        # We cover all the vertices in the MTG, from the tips to the base:
        for vid in post_order(self.g, root):

            # n represents the current root element:
            n = self.g.node(vid)

            # We calculate a coefficient that will modify the different "ages" experienced by roots according to soil
            # temperature assuming a linear relationship (this is equivalent as the calculation of "growth degree-days):
            temperature_time_adjustment = self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=n.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)
    
            # AVOIDANCE OF UNWANTED CASES:
            # -----------------------------
            # We make sure that the element is not dead:
            if n.type == "Dead" or n.type == "Just_dead" or n.type == "Support_for_seminal_root" or n.type == "Support_for_adventitious_root":
                # In such case, we just pass to the next element in the iteration:
                continue

            # We make sure that there is a potential growth for this element:
            if n.potential_length <= n.initial_length and n.potential_radius <= n.initial_radius:
                # In such case, we just pass to the next element in the iteration:
                continue
            # INITIALIZATION AND CALCULATIONS OF POTENTIAL GROWTH DEMAND IN HEXOSE:
            # ----------------------------------------------------------------------

            # WARNING: All growth related variables should have been initialized by another module at the beginning of the time step!!!

            # We calculate the initial volume of the element:
            initial_volume = self.volume_from_radius_and_length(n, n.initial_radius, n.initial_length)
            # We calculate the potential volume of the element based on the potential radius and potential length:
            potential_volume = self.volume_from_radius_and_length(n, n.potential_radius, n.potential_length)
            # We calculate the number of moles of hexose required for growth, including the respiration cost according to
            # the yield growth included in the model of Thornley and Cannell (2000), where root_tissue_density is the dry structural
            # weight per volume (g m-3) and struct_mass_C_content is the amount of C per gram of dry structural mass (mol_C g-1):
            n.hexose_growth_demand = (potential_volume - initial_volume) \
                                     * n.root_tissue_density * self.struct_mass_C_content / self.yield_growth * 1 / 6.
            # We verify that this potential growth demand is positive:
            if n.hexose_growth_demand < 0.:
                print("!!! ERROR: a negative growth demand of", n.hexose_growth_demand,
                      "was calculated for the element", n.index(), "of class", n.label)
                print("The initial volume is", initial_volume, "the potential volume is", potential_volume)
                print("The initial length was", n.initial_length, "and the potential length was",
                      n.potential_length)
                print("The initial radius was", n.initial_radius, "and the potential radius was",
                      n.potential_radius)
                n.hexose_growth_demand = 0.
                # In such case, we just pass to the next element in the iteration:
                continue
            elif n.hexose_growth_demand == 0.:
                continue

            n.amino_acids_growth_demand = (potential_volume - initial_volume) \
                                     * n.root_tissue_density * max(self.struct_mass_N_content / self.yield_growth_N / self.r_Nm_AA,
                                                                # We also increase the demand as amino acids also bring C!
                                                                self.struct_mass_C_content / self.yield_growth / self.r_C_AA)
            # We verify that this potential growth demand is positive:
            if n.amino_acids_growth_demand < 0.:
                print("!!! ERROR: a negative growth demand for amino acids of", n.amino_acids_growth_demand,
                      "was calculated for the element", n.index(), "of class", n.label)
                print("The initial volume is", initial_volume, "the potential volume is", potential_volume)
                print("The initial length was", n.initial_length, "and the potential length was",
                      n.potential_length)
                print("The initial radius was", n.initial_radius, "and the potential radius was",
                      n.potential_radius)
                n.amino_acids_growth_demand = 0.
                # In such case, we just pass to the next element in the iteration:
                continue
            elif n.amino_acids_growth_demand == 0.:
                continue

            # CALCULATIONS OF THE AMOUNT OF HEXOSE AVAILABLE FOR GROWTH:
            # ---------------------------------------------------------

            # We initialize each amount of hexose available for growth:
            hexose_available_for_elongation = 0.
            hexose_available_for_thickening = 0.
            amino_acids_possibly_required_for_elongation = 0.
            amino_acids_available_for_thickening = 0.

            # If elongation is possible:
            if n.potential_length > n.length:
                # TODO: Check how useful it is to limit the growth so that the concentration does not go below a treshold.
                # NEW - We calculate the actual amount of hexose that can be used for growth according to the treshold
                # concentration below which no growth process is authorized:
                hexose_available_for_elongation \
                    = n.hexose_possibly_required_for_elongation \
                    - n.struct_mass_contributing_to_elongation * self.C_hexose_min_for_elongation
                if hexose_available_for_elongation < 0.:
                    hexose_available_for_elongation = 0.
                amino_acids_possibly_required_for_elongation = n.amino_acids_possibly_required_for_elongation
                list_of_elongation_supporting_elements = n.list_of_elongation_supporting_elements
                list_of_elongation_supporting_elements_hexose = n.list_of_elongation_supporting_elements_hexose
                list_of_elongation_supporting_elements_amino_acids = n.list_of_elongation_supporting_elements_amino_acids
                list_of_elongation_supporting_elements_mass = n.list_of_elongation_supporting_elements_mass
            
            # If radial growth is possible:
            if n.potential_radius > n.radius:
                # We only consider the amount of hexose immediately available in the element that can increase in radius:
                # We only consider the amount of hexose immediately available in the element that can increase in radius.
                # TODO: Check how useful it is to limit the growth so that the concentration does not go below a treshold.
                # NEW - We calculate the actual amount of hexose that can be used for thickening according to the treshold
                # concentration below which no thickening is authorized:
                hexose_available_for_thickening = n.hexose_available_for_thickening \
                                                - n.struct_mass * self.C_hexose_min_for_thickening
                if hexose_available_for_thickening <0.:
                    hexose_available_for_thickening = 0.

                amino_acids_available_for_thickening = n.amino_acids_available_for_thickening

            # In case no hexose is available at all:
            if (hexose_available_for_elongation  + hexose_available_for_thickening) <= 0. or (
                amino_acids_possibly_required_for_elongation + amino_acids_available_for_thickening) <= 0. :
                # Then we move to the next element in the main loop:
                continue

            # We initialize the temporary variable "remaining_hexose" that computes the amount of hexose left for growth:
            remaining_hexose_for_elongation = hexose_available_for_elongation 
            remaining_hexose_for_thickening = hexose_available_for_thickening
            remaining_amino_acids_for_elongation = amino_acids_possibly_required_for_elongation
            remaining_amino_acids_for_thickening = amino_acids_available_for_thickening

            # ACTUAL ELONGATION IS FIRST CONSIDERED:
            # ---------------------------------------

            # We calculate the maximal possible length of the root element according to all the hexose available for elongation:
            prossible_volume_increase_C_hexose = hexose_available_for_elongation * 6. \
                         / (n.root_tissue_density * self.struct_mass_C_content) * self.yield_growth
            prossible_volume_increase_C_AA = amino_acids_possibly_required_for_elongation * self.r_C_AA \
                         / (n.root_tissue_density * self.struct_mass_C_content) * self.yield_growth
            volume_max_C = initial_volume + prossible_volume_increase_C_hexose + prossible_volume_increase_C_AA
            volume_max_N = initial_volume + amino_acids_possibly_required_for_elongation * self.r_Nm_AA \
                         / (n.root_tissue_density * self.struct_mass_N_content) * self.yield_growth_N
            # We account for the minimal volume defining the most limiting factor between C and N
            length_max = min(volume_max_C, volume_max_N) / (pi * n.initial_radius ** 2)

            # If the element can elongate:
            if n.potential_length > n.initial_length:
                # CALCULATING ACTUAL ELONGATION:
                # If elongation is possible but is limited by the amount of hexose available:
                if n.potential_length >= length_max:
                    # Elongation is limited using all the amount of hexose available:
                    n.length = length_max
                # Otherwise, elongation can be done up to the full potential:
                else:
                    # Elongation is done up to the full potential:
                    n.length = n.potential_length
                # The corresponding new volume is calculated:
                volume_after_elongation = self.volume_from_radius_and_length(n, n.initial_radius, n.length)

                # The overall cost of elongation is calculated as:
                amino_acids_consumption_by_elongation = (1. / self.r_Nm_AA) * (volume_after_elongation - initial_volume) \
                    * n.root_tissue_density * self.struct_mass_N_content / self.yield_growth_N

                # Here we substract the amount of C brought by amino acids to the amount of consummed hexose
                hexose_consumption_by_elongation = (1. / 6. * (volume_after_elongation - initial_volume) \
                    * n.root_tissue_density * self.struct_mass_C_content / self.yield_growth) - amino_acids_consumption_by_elongation * self.r_C_AA / 6

                # If there has been an actual elongation:
                if n.length > n.initial_length:

                    # REGISTERING THE COSTS FOR ELONGATION:
                    # We cover each of the elements that have provided hexose for sustaining the elongation of element n:
                    for i in range(0, len(list_of_elongation_supporting_elements)):
                        index = list_of_elongation_supporting_elements[i]
                        supplying_element = self.g.node(index)
                        # We define the actual contribution of the current element based on total hexose consumption by growth
                        # of element n and the relative contribution of the current element to the pool of the potentially available hexose:
                        hexose_actual_contribution_to_elongation = hexose_consumption_by_elongation \
                                                                   * list_of_elongation_supporting_elements_hexose[
                                                                       i] / n.hexose_possibly_required_for_elongation
                        amino_acids_actual_contribution_to_elongation = amino_acids_consumption_by_elongation \
                                                                   * list_of_elongation_supporting_elements_amino_acids[
                                                                       i] / amino_acids_possibly_required_for_elongation
                        # The amount of hexose used for growth in this element is increased:
                        supplying_element.hexose_consumption_by_growth_amount += hexose_actual_contribution_to_elongation
                        supplying_element.hexose_consumption_by_growth += hexose_actual_contribution_to_elongation / self.time_step_in_seconds
                        # The amount of amino acids used for growth in this element is increased:
                        supplying_element.amino_acids_consumption_by_growth_amount += amino_acids_actual_contribution_to_elongation
                        supplying_element.amino_acids_consumption_by_growth += amino_acids_actual_contribution_to_elongation / self.time_step_in_seconds
                        # And the amount of hexose that has been used for growth respiration is calculated and transformed into moles of CO2:
                        supplying_element.resp_growth += hexose_actual_contribution_to_elongation * (1 - self.yield_growth) * 6.

            # ACTUAL RADIAL GROWTH IS THEN CONSIDERED:
            # -----------------------------------------
            # If the radius of the element can increase:
            if n.potential_radius > n.initial_radius:
                # CALCULATING ACTUAL THICKENING:
                # We calculate the increase in volume that can be achieved with the amount of hexose available:
                possible_radial_increase_in_volume_C_hexose = \
                    remaining_hexose_for_thickening * 6. * self.yield_growth \
                    / (n.root_tissue_density * self.struct_mass_C_content)
                possible_radial_increase_in_volume_C_AA = \
                    remaining_amino_acids_for_thickening * self.r_C_AA * self.yield_growth \
                    / (n.root_tissue_density * self.struct_mass_C_content)
                possible_radial_increase_in_volume_C = possible_radial_increase_in_volume_C_hexose + possible_radial_increase_in_volume_C_AA
                possible_radial_increase_in_volume_N = \
                    remaining_amino_acids_for_thickening * self.r_Nm_AA * self.yield_growth_N \
                    / (n.root_tissue_density * self.struct_mass_N_content)
                # We calculate the maximal possible volume based on the volume of the new cylinder after elongation
                # and the increase in volume that could be achieved by consuming the first limiting factor between hexose and amino acids:

                volume_max = self.volume_from_radius_and_length(n, n.initial_radius, n.length) + min(possible_radial_increase_in_volume_C, possible_radial_increase_in_volume_N)
                # We then calculate the corresponding new possible radius corresponding to this maximum volume:
                if n.type == "Root_nodule":
                    # If the element corresponds to a nodule, then it we calculate the radius of a theoretical sphere:
                    possible_radius = (3. / (4. * pi)) ** (1. / 3.)
                else:
                    # Otherwise, we calculate the radius of a cylinder:
                    possible_radius = sqrt(volume_max / (n.length * pi))
                if possible_radius < 0.9999 * n.initial_radius:  # We authorize a difference of 0.01% due to calculation errors!
                    print("!!! ERROR: the calculated new radius of element", n.index(),
                          "is lower than the initial one!")
                    print("The possible radius was", possible_radius, "and the initial radius was",
                          n.initial_radius)

                # If the maximal radius that can be obtained is lower than the potential radius suggested by the potential growth module:
                if possible_radius <= n.potential_radius:
                    # Then radial growth is limited and there is no remaining hexose after radial growth:
                    n.radius = possible_radius
                    hexose_actual_contribution_to_thickening = remaining_hexose_for_thickening
                    remaining_hexose_for_thickening = 0.
                    amino_acids_actual_contribution_to_thickening = remaining_amino_acids_for_thickening
                    remaining_amino_acids_for_thickening = 0.
                else:
                    # Otherwise, radial growth is done up to the full potential and the remaining hexose is calculated:
                    n.radius = n.potential_radius
                    net_increase_in_volume = self.volume_from_radius_and_length(n, n.radius, n.length) \
                        - self.volume_from_radius_and_length(n, n.initial_radius, n.length)
                    # net_increase_in_volume = pi * (n.radius ** 2 - n.initial_radius ** 2) * n.length
                    # We then calculate the remaining amount of hexose after thickening:

                    amino_acids_actual_contribution_to_thickening = \
                        1. / self.r_Nm_AA * net_increase_in_volume \
                        * n.root_tissue_density * self.struct_mass_N_content / self.yield_growth_N

                    hexose_actual_contribution_to_thickening = \
                        1. / 6. * net_increase_in_volume \
                        * n.root_tissue_density * self.struct_mass_C_content / self.yield_growth - amino_acids_actual_contribution_to_thickening * self.r_C_AA / 6
                    

                # REGISTERING THE COSTS FOR THICKENING:
                # --------------------------------------
                fraction_of_available_hexose_in_the_element = \
                    (n.C_hexose_root * n.initial_struct_mass) / hexose_available_for_thickening
                # The amount of hexose used for growth in this element is increased:
                n.hexose_consumption_by_growth_amount += \
                    (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element)
                n.hexose_consumption_by_growth += \
                    (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element) / self.time_step_in_seconds
                # Same calculation for amino acids costs
                fraction_of_available_amino_acids_in_the_element = \
                    (n.AA * n.initial_struct_mass) / amino_acids_available_for_thickening
                # The amount of amino acids used for growth in this element is increased:
                n.amino_acids_consumption_by_growth_amount += \
                    (amino_acids_actual_contribution_to_thickening * fraction_of_available_amino_acids_in_the_element)
                n.hexose_consumption_by_growth += \
                    (amino_acids_actual_contribution_to_thickening * fraction_of_available_amino_acids_in_the_element) / self.time_step_in_seconds
                # And the amount of hexose that has been used for growth respiration is calculated and transformed into moles of CO2:
                n.resp_growth += \
                    (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element) \
                    * (1 - self.yield_growth) * 6.
                if n.type == "Root_nodule":
                    index_parent = self.g.Father(n.index(), EdgeType='+')
                    parent = self.g.node(index_parent)
                    fraction_of_available_hexose_in_the_element = \
                        (parent.C_hexose_root * parent.initial_struct_mass) / hexose_available_for_thickening
                    # The amount of hexose used for growth in this element is increased:
                    parent.hexose_consumption_by_growth_amount += \
                        (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element)
                    parent.hexose_consumption_by_growth += \
                        (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element) / self.time_step_in_seconds
                    # Same calculation for amino acids costs
                    fraction_of_available_amino_acids_in_the_element = \
                        (parent.AA * parent.initial_struct_mass) / amino_acids_available_for_thickening
                    # The amount of hexose used for growth in this element is increased:
                    parent.amino_acids_consumption_by_growth_amount += \
                        (amino_acids_actual_contribution_to_thickening * fraction_of_available_amino_acids_in_the_element)
                    parent.amino_acids_consumption_by_growth += \
                        (amino_acids_actual_contribution_to_thickening * fraction_of_available_amino_acids_in_the_element) / self.time_step_in_seconds
                    
                    # And the amount of hexose that has been used for growth respiration is calculated and transformed into moles of CO2:
                    parent.resp_growth += \
                        (hexose_actual_contribution_to_thickening * fraction_of_available_hexose_in_the_element) \
                        * (1 - self.yield_growth) * 6.

            # RECORDING THE ACTUAL STRUCTURAL MODIFICATIONS:
            # -----------------------------------------------
            # The new volume of the element is automatically calculated
            n.volume = self.volume_from_radius_and_length(n, n.radius, n.length)
            # The new dry structural struct_mass of the element is calculated from its new volume:
            n.struct_mass = n.volume * n.root_tissue_density
            n.struct_mass_produced = (n.volume - initial_volume) * n.root_tissue_density

            if n.struct_mass < n.initial_struct_mass and n.struct_mass_produced > 0.:
                print(f"!!! ERROR during initialisation for initial struct mass, no concentrations will be updated on {n.index()}")
                n.initial_struct_mass = n.struct_mass

            # Verification: we check that no negative length or struct_mass have been generated!
            if n.volume < 0:
                print("!!! ERROR: the element", n.index(), "of class", n.label, "has a length of", n.length,
                      "and a mass of", n.struct_mass)
                # We then reset all the geometrical values to their initial values:
                n.length = n.initial_length
                n.radius = n.initial_radius
                n.struct_mass = n.initial_struct_mass
                n.struct_mass_produced = 0.
                n.volume = initial_volume

            # MODIFYING SPECIFIC PROPERTIES AFTER ELONGATION:
            # -----------------------------------------------
            # If there has been an actual elongation:
            if n.length > n.initial_length:
                # If the elongated apex corresponded to any primordium that has been allowed to emerge:
                if n.type == "Seminal_root_before_emergence" \
                        or n.type == "Adventitious_root_before_emergence" \
                        or n.type == "Normal_root_before_emergence":
                    # We now consider the apex to have emerged:
                    n.type = "Normal_root_after_emergence"
                    # The exact time since emergence is recorded:
                    n.thermal_time_since_emergence = n.thermal_potential_time_since_emergence
                    n.actual_time_since_emergence = n.thermal_time_since_emergence / temperature_time_adjustment
                    # The actual elongation rate is calculated:
                    n.actual_elongation = n.length - n.initial_length
                    n.actual_elongation_rate = n.actual_elongation / n.actual_time_since_emergence
                    # Note: at this stage, no sugar has been allocated to the emerging primordium itself!
                    # if n.type == "Adventitious_root_before_emergence":
                    #     print("> A new adventitious root has emerged, starting from element", n.index(), "!")
                elif n.type == "Normal_root_after_emergence":
                    # The actual elongation rate is calculated:
                    n.actual_elongation = n.length - n.initial_length
                    n.actual_elongation_rate = n.actual_elongation / self.time_step_in_seconds

                # # IF THE AGE OF CELLS NEEDS TO BE UPDATED AT THIS STAGE:
                # # In both cases, we calculate a new, average age of the root cells in the elongated apex, based on the
                # # principle that cells appear at the very end of the root tip, and then age. In this case, the average age
                # # of the element that has elongated is calculated as the mean value between the incremented age of the part
                # # that was already formed and the average age of the new part that has been created:
                # initial_time_since_cells_formation = n.actual_time_since_cells_formation
                # Age_non_elongated = initial_time_since_cells_formation + time_step_in_seconds
                # Age_elongated =  0.5 * n.actual_elongation / n.actual_elongation_rate
                # n.actual_time_since_cells_formation = (n.initial_length * Age_non_elongated + n.actual_elongation * Age_elongated) / n.length
                # n.thermal_time_since_cells_formation += (n.actual_time_since_cells_formation - initial_time_since_cells_formation) * temperature_time_adjustment

                # The distance to the last ramification is increased:
                n.dist_to_ramif += n.actual_elongation


    @postsegmentation
    @state
    def root_hairs_dynamics(self):
        
        """
        This function computes the evolution of the density and average length of root hairs along each root,
        and specifies which hairs are alive or dead.

        EDIT : Added a limitation of growth for root hairs associated to the availability of amino acids availability in the root.
        Accounted for the consumption in the consumption of amino acids by growth.
        :return:
        """

        #  TODO FOR TRISTAN: In a second step, consider playing on the density / max. length of root hairs depending on the availability of N in the soil (if relevant)?

        # We cover all the vertices in the MTG:
        for vid in self.g.vertices_iter(scale=1):
            # n represents the vertex:
            n = self.g.node(vid)

            # First, we ensure that the element has a positive length:
            if n.length <= 0:
                continue

            # We also exclude nodules and dead elements from this computation:
            if n.type == "Just_dead" or n.type == "Dead" or n.type == "Nodule":
                continue

            # # TODO: Check the consequences of avoiding apex in root hairs dynamics!
            # # WE ALSO AVOID ROOT APICES - EVEN IF IN THEORY ROOT HAIRS MAY ALSO APPEAR ON THEM:
            # if n.label == "Apex":
            #     continue
            # # Even if root hairs should have already emerge on that root apex, they will appear in the next step (or in a few steps)
            # # when the element becomes a segment.

            # We calculate the equivalent of a thermal time for the current time step:
            temperature_time_adjustment = self.temperature_modification(process_at_T_ref=self.process_at_T_ref,
                                                                    soil_temperature=n.soil_temperature,
                                                                    T_ref=self.T_ref, A=self.A, B=self.B, C=self.C)
            elapsed_thermal_time = self.time_step_in_seconds * temperature_time_adjustment

            # We keep in memory the initial total mass of root hairs (possibly including dead hairs):
            initial_root_hairs_struct_mass = n.root_hairs_struct_mass

            # We calculate the total number of (newly formed) root hairs (if any) and update their age:
            # ------------------------------------------------------------------------------------------
            # CASE 1 - If the current element is completely included within the actual growing zone of the root at the root
            # tip, the root hairs cannot have formed yet:
            if n.distance_from_tip <= self.growing_zone_factor * n.radius:
                # We stop here with the calculations and move to the next element:
                continue
            # CASE 2 - If all root hairs have already been formed:
            if n.all_root_hairs_formed:
                # Then we simply increase the time since root hairs emergence started:
                n.actual_time_since_root_hairs_emergence_started += self.time_step_in_seconds
                n.thermal_time_since_root_hairs_emergence_started += elapsed_thermal_time
                n.actual_time_since_root_hairs_emergence_stopped += self.time_step_in_seconds
                n.thermal_time_since_root_hairs_emergence_stopped += elapsed_thermal_time
            # CASE 3 - If the theoretical growing zone limit is located somewhere within the root element:
            elif n.distance_from_tip - n.length < self.growing_zone_factor * n.radius:
                # We first record the previous length of the root hair zone within the element:
                initial_length_with_hairs = n.actual_length_with_hairs
                # Then the new length of the root hair zone is calculated:
                n.actual_length_with_hairs = n.distance_from_tip - self.growing_zone_factor * n.radius
                net_increase_in_root_hairs_length = n.actual_length_with_hairs - initial_length_with_hairs
                # The corresponding number of root hairs is calculated:
                n.total_root_hairs_number = self.root_hairs_density * n.radius * n.actual_length_with_hairs
                # The time since root hair formation started is then calculated, using the recent increase in the length
                # of the current root hair zone and the elongation rate of the corresponding root tip. The latter is
                # calculated using the difference between the new distance_from_tip of the element and the previous one:
                elongation_rate_in_actual_time = (n.distance_from_tip - n.former_distance_from_tip) / self.time_step_in_seconds
                elongation_rate_in_thermal_time = (n.distance_from_tip - n.former_distance_from_tip) / elapsed_thermal_time
                # SUBCASE 3.1 - If root hairs had not emerged at the previous time step:
                if elongation_rate_in_actual_time > 0. and initial_length_with_hairs <= 0.:
                    # We increase the time since root hairs emerged by only the fraction of the time step corresponding to the growth of hairs:
                    n.actual_time_since_root_hairs_emergence_started += \
                        self.time_step_in_seconds - net_increase_in_root_hairs_length / elongation_rate_in_actual_time
                    n.thermal_time_since_root_hairs_emergence_started += \
                        elapsed_thermal_time - net_increase_in_root_hairs_length / elongation_rate_in_thermal_time
                # SUBCASE 3.2 - the hairs had already started to grow:
                else:
                    # Consequently, the full time elapsed during this time step can be added to the age:
                    n.actual_time_since_root_hairs_emergence_started += self.time_step_in_seconds
                    n.thermal_time_since_root_hairs_emergence_started += elapsed_thermal_time
            # CASE 4 - the element is now "full" with root hairs as the limit of root elongation is located further down:
            else:
                # The actual time since root hairs emergence started is first increased:
                n.actual_time_since_root_hairs_emergence_started += self.time_step_in_seconds
                n.thermal_time_since_root_hairs_emergence_started += elapsed_thermal_time
                # We then record the previous length of the root hair zone within the root element:
                initial_length_with_hairs = n.actual_length_with_hairs
                # And the new length of the root hair zone is necessarily the full length of the root element:
                n.actual_length_with_hairs = n.length
                net_increase_in_root_hairs_length = n.actual_length_with_hairs - initial_length_with_hairs
                # The total number of hairs is defined according to the radius and total length of the element:
                n.total_root_hairs_number = self.root_hairs_density * n.radius * n.length
                # The elongation of the corresponding root tip is calculated as the difference between the new
                # distance_from_tip of the element and the previous one:
                elongation_rate_in_actual_time = (
                                                         n.distance_from_tip - n.former_distance_from_tip) / self.time_step_in_seconds
                elongation_rate_in_thermal_time = (
                                                          n.distance_from_tip - n.former_distance_from_tip) / elapsed_thermal_time
                # The actual time since root hairs emergence has stopped is then calculated:
                if elongation_rate_in_actual_time > 0.:
                    n.actual_time_since_root_hairs_emergence_stopped += \
                        self.time_step_in_seconds - net_increase_in_root_hairs_length / elongation_rate_in_actual_time
                    n.thermal_time_since_root_hairs_emergence_stopped += \
                        elapsed_thermal_time - net_increase_in_root_hairs_length / elongation_rate_in_thermal_time
                else:
                    n.actual_time_since_root_hairs_emergence_stopped += self.time_step_in_seconds
                    n.thermal_time_since_root_hairs_emergence_stopped += elapsed_thermal_time
                # At this stage, all root hairs that could be formed have been formed, so we record this:
                n.all_root_hairs_formed = True

            # We now calculate the number of living and dead root hairs:
            # -----------------------------------------------------------
            # Root hairs are dying when the time since they emerged is higher than their lifespan.
            # If the time since root hairs emergence started is lower than the lifespan,
            # no root hair should be dead:
            if n.thermal_time_since_root_hairs_emergence_started <= n.root_hairs_lifespan:
                n.dead_root_hairs_number = 0.
            # Otherwise, if the time since root hairs emergence stopped is higher than the lifespan:
            elif n.thermal_time_since_root_hairs_emergence_stopped > n.root_hairs_lifespan:
                # Then all the root hairs of the root element must now be dead:
                n.dead_root_hairs_number = n.total_root_hairs_number
            # In the intermediate case, there are currently both dead and living root hairs on the root element:
            else:
                # We assume that there is a linear decrease of root hair age between the first hair that has emerged
                # and the last one that has emerged:
                time_since_first_death = n.thermal_time_since_root_hairs_emergence_started - n.root_hairs_lifespan
                dead_fraction = time_since_first_death / (n.thermal_time_since_root_hairs_emergence_started
                                                          - n.thermal_time_since_root_hairs_emergence_stopped)
                n.dead_root_hairs_number = n.total_root_hairs_number * dead_fraction

            # In all cases, the number of the living root hairs is then calculated by difference with the total hair number:
            n.living_root_hairs_number = n.total_root_hairs_number - n.dead_root_hairs_number

            # We calculate the new average root hairs length, if needed:
            # ----------------------------------------------------------
            # If the root hairs had not reached their maximal length:
            if n.root_hair_length < self.root_hair_max_length:
                # The new potential root hairs length is calculated according to the elongation rate,
                # corrected by temperature and modulated by the concentration of hexose (in the same way as for root
                # elongation) available in the root hair zone on the root element:
                cn_limitation = ((n.C_hexose_root / (n.C_hexose_root + self.Km_elongation)) + (n.AA / (n.AA + self.Km_elongation_amino_acids))) / 2
                new_length = n.root_hair_length + self.root_hairs_elongation_rate * self.root_hair_radius * (n.actual_length_with_hairs / n.length) * cn_limitation * elapsed_thermal_time
                
                # If the new calculated length is higher than the maximal length:
                if new_length > self.root_hair_max_length:
                    # We set the root hairs length to the maximal length:
                    n.root_hair_length = self.root_hair_max_length
                else:
                    # Otherwise, we record the new calculated length:
                    n.root_hair_length = new_length

            # We finally calculate the total external surface (m2), volume (m3) and mass (g) of root hairs:
            # ----------------------------------------------------------------------------------------------
            # In the calculation of surface, we consider the root hair to be a cylinder, and include the lateral section,
            # but exclude the section of the cylinder at the tip:
            n.root_hairs_volume = (self.root_hair_radius ** 2 * pi) * n.root_hair_length * n.total_root_hairs_number
            n.root_hairs_struct_mass = n.root_hairs_volume * n.root_tissue_density
            if n.total_root_hairs_number > 0.:
                n.living_root_hairs_struct_mass = n.root_hairs_struct_mass * n.living_root_hairs_number \
                                                  / n.total_root_hairs_number
            else:
                n.living_root_hairs_struct_mass = 0.

            # We calculate the mass of hairs that has been effectively produced, including from root hairs that may have died since then:
            # ----------------------------------------------------------------------------------------------------------------------------
            # We calculate the new production as the difference between initial and final mass:
            n.root_hairs_struct_mass_produced = n.root_hairs_struct_mass - initial_root_hairs_struct_mass

            # We add the cost of producing the new living root hairs (if any) to the hexose consumption by growth:
            hexose_consumption = n.root_hairs_struct_mass_produced * self.struct_mass_C_content / self.yield_growth / 6.
            amino_acids_comsumption = n.root_hairs_struct_mass_produced * self.struct_mass_N_content / self.yield_growth_N / self.r_Nm_AA

            n.hexose_consumption_by_growth_amount += hexose_consumption
            n.hexose_consumption_by_growth += hexose_consumption / self.time_step_in_seconds
            n.amino_acids_consumption_by_growth_amount += amino_acids_comsumption
            n.amino_acids_consumption_by_growth += amino_acids_comsumption / self.time_step_in_seconds
            n.resp_growth += hexose_consumption * 6. * (1 - self.yield_growth)

    # Adding a new root element with pre-defined properties:
    def ADDING_A_CHILD(self, mother_element, edge_type='+', label='Apex', type='Normal_root_before_emergence',
                       root_order=1, angle_down=45., angle_roll=0., length=0., radius=0.,
                       identical_properties=True, nil_properties=False):
        """
        This function creates a new child element on the mother element, based on the function add_child.
        When called, this allows to automatically define standard properties without defining them in the main code.
        :param mother_element: the node of the MTG on which the child element will be created
        :param edge_type: the type of relationship between the mother element and its child ('+' or '<')
        :param label: label of the child element
        :param type: type of the child element
        :param root_order: root_order of the child element
        :param angle_down: angle_down of the child element
        :param angle_roll: angle_roll of the child element
        :param length: length of the child element
        :param radius: radius of the child element
        :param identical_properties: if True, the main properties of the child will be identical to those of the mother
        :param nil_properties: if True, the main properties of the child will be 0
        :return: the new child element
        """

        # TODO# FOR TRISTAN: When working with a dynamic root structure, you will need to specify in this function
        #  "ADDING_A_CHILD" your new variables that will either be set to 0 (nil properties) or be equal to that of the mother
        #  element.

        # If nil_properties = True, then we set most of the properties of the new element to 0:
        if nil_properties:
            new_child = mother_element.add_child(edge_type=edge_type,
                                                 # Characteristics:
                                                 # -----------------
                                                 label=label,
                                                 type=type,
                                                 root_order=root_order,
                                                 # Authorizations and C requirements:
                                                 # -----------------------------------
                                                 lateral_root_emergence_possibility='Impossible',
                                                 emergence_cost=0.,
                                                 # Geometry and topology:
                                                 # -----------------------
                                                 angle_down=angle_down,
                                                 angle_roll=angle_roll,
                                                 # The length of the primordium is set to 0:
                                                 length=length,
                                                 radius=radius,
                                                 original_radius=radius,
                                                 potential_length=0.,
                                                 theoretical_radius=radius,
                                                 potential_radius=radius,
                                                 initial_length=0.,
                                                 initial_radius=radius,
                                                 volume=0.,
                                                 root_tissue_density=self.new_root_tissue_density,
                                                 dist_to_ramif=0.,
                                                 distance_from_tip=0.,
                                                 former_distance_from_tip=0.,
                                                 actual_elongation=0.,
                                                 actual_elongation_rate=0.,
                                                 # Quantities and concentrations:
                                                 # -------------------------------
                                                 struct_mass=0.,
                                                 initial_struct_mass=0.,
                                                 C_hexose_root=0.,
                                                 # Root hairs:
                                                 # ------------
                                                 root_hair_radius=self.root_hair_radius,
                                                 root_hair_length=0.,
                                                 actual_length_with_hairs=0.,
                                                 living_root_hairs_number=0.,
                                                 dead_root_hairs_number=0.,
                                                 total_root_hairs_number=0.,
                                                 actual_time_since_root_hairs_emergence_started=0.,
                                                 thermal_time_since_root_hairs_emergence_started=0.,
                                                 actual_time_since_root_hairs_emergence_stopped=0.,
                                                 thermal_time_since_root_hairs_emergence_stopped=0.,
                                                 all_root_hairs_formed=False,
                                                 root_hairs_lifespan=self.root_hairs_lifespan,
                                                 root_hairs_struct_mass=0.,
                                                 root_hairs_struct_mass_produced=0.,
                                                 initial_living_root_hairs_struct_mass=0.,
                                                 living_root_hairs_struct_mass=0.,
                                                 # Fluxes:
                                                 # --------
                                                 resp_growth=0.,
                                                 struct_mass_produced=0.,
                                                 hexose_growth_demand=0.,
                                                 hexose_consumption_by_growth_amount=0.,
                                                 hexose_consumption_by_growth=0.,
                                                 hexose_consumption_by_fungus=0.,
                                                 hexose_possibly_required_for_elongation=0.,
                                                 # Time indications:
                                                 # ------------------
                                                 soil_temperature=7.8,  # TODO change
                                                 growth_duration=self.GDs * (2 * radius) ** 2,
                                                 life_duration=self.LDs * 2. * radius * self.new_root_tissue_density,
                                                 actual_time_since_primordium_formation=0.,
                                                 actual_time_since_emergence=0.,
                                                 actual_time_since_cells_formation=0.,
                                                 actual_potential_time_since_emergence=0.,
                                                 actual_time_since_growth_stopped=0.,
                                                 actual_time_since_death=0.,
                                                 thermal_time_since_primordium_formation=0.,
                                                 thermal_time_since_emergence=0.,
                                                 thermal_time_since_cells_formation=0.,
                                                 thermal_potential_time_since_emergence=0.,
                                                 thermal_time_since_growth_stopped=0.,
                                                 thermal_time_since_death=0.,
                                                 AA=0.,
                                                 amino_acids_consumption_by_growth_amount=0.,
                                                 amino_acids_consumption_by_growth=0.,
                                                 amino_acids_possibly_required_for_elongation=0.,
                                                 amino_acids_growth_demand=0.
                                                 )
            
            return new_child

        # Otherwise, if identical_properties=True, then we copy most of the properties of the mother element in the new element:
        elif identical_properties:
            new_child = mother_element.add_child(edge_type=edge_type,
                                                 # Characteristics:
                                                 # -----------------
                                                 label=label,
                                                 type=type,
                                                 root_order=root_order,
                                                 # Authorizations and C requirements:
                                                 # -----------------------------------
                                                 lateral_root_emergence_possibility='Impossible',
                                                 emergence_cost=0.,
                                                 # Geometry and topology:
                                                 # -----------------------
                                                 angle_down=angle_down,
                                                 angle_roll=angle_roll,
                                                 # The length of the primordium is set to 0:
                                                 length=length,
                                                 radius=mother_element.radius,
                                                 original_radius=mother_element.radius,
                                                 potential_length=length,
                                                 theoretical_radius=mother_element.theoretical_radius,
                                                 potential_radius=mother_element.potential_radius,
                                                 initial_length=length,
                                                 initial_radius=mother_element.radius,
                                                 volume=0.,
                                                 root_tissue_density=mother_element.root_tissue_density,
                                                 dist_to_ramif=mother_element.dist_to_ramif,
                                                 distance_from_tip=mother_element.distance_from_tip,
                                                 former_distance_from_tip=mother_element.former_distance_from_tip,
                                                 actual_elongation=mother_element.actual_elongation,
                                                 actual_elongation_rate=mother_element.actual_elongation_rate,
                                                 # Quantities and concentrations:
                                                 # -------------------------------
                                                 struct_mass=mother_element.struct_mass,
                                                 initial_struct_mass=mother_element.initial_struct_mass,
                                                 C_hexose_root=mother_element.C_hexose_root,
                                                 # Root hairs:
                                                 # ------------
                                                 root_hair_radius=mother_element.root_hair_radius,
                                                 root_hair_length=mother_element.root_hair_length,
                                                 actual_length_with_hairs=mother_element.actual_length_with_hairs,
                                                 living_root_hairs_number=mother_element.living_root_hairs_number,
                                                 dead_root_hairs_number=mother_element.dead_root_hairs_number,
                                                 total_root_hairs_number=mother_element.total_root_hairs_number,
                                                 actual_time_since_root_hairs_emergence_started=mother_element.actual_time_since_root_hairs_emergence_started,
                                                 thermal_time_since_root_hairs_emergence_started=mother_element.thermal_time_since_root_hairs_emergence_started,
                                                 actual_time_since_root_hairs_emergence_stopped=mother_element.actual_time_since_root_hairs_emergence_stopped,
                                                 thermal_time_since_root_hairs_emergence_stopped=mother_element.thermal_time_since_root_hairs_emergence_stopped,
                                                 all_root_hairs_formed=mother_element.all_root_hairs_formed,
                                                 root_hairs_lifespan=mother_element.root_hairs_lifespan,
                                                 root_hairs_struct_mass=mother_element.root_hairs_struct_mass,
                                                 root_hairs_struct_mass_produced=mother_element.root_hairs_struct_mass_produced,
                                                 living_root_hairs_struct_mass=mother_element.living_root_hairs_struct_mass,
                                                 initial_living_root_hairs_struct_mass=mother_element.initial_living_root_hairs_struct_mass,
                                                 # Fluxes:
                                                 # -------
                                                 resp_growth=mother_element.resp_growth,
                                                 struct_mass_produced=mother_element.struct_mass_produced,
                                                 hexose_growth_demand=mother_element.hexose_growth_demand,
                                                 hexose_possibly_required_for_elongation=mother_element.hexose_possibly_required_for_elongation,
                                                 hexose_consumption_by_growth_amount=mother_element.hexose_consumption_by_growth_amount,
                                                 hexose_consumption_by_growth=mother_element.hexose_consumption_by_growth,
                                                 hexose_consumption_by_fungus=mother_element.hexose_consumption_by_fungus,
                                                 # Time indications:
                                                 # ------------------
                                                 soil_temperature=mother_element.soil_temperature,
                                                 growth_duration=mother_element.growth_duration,
                                                 life_duration=mother_element.life_duration,
                                                 actual_time_since_primordium_formation=mother_element.actual_time_since_primordium_formation,
                                                 actual_time_since_emergence=mother_element.actual_time_since_emergence,
                                                 actual_time_since_cells_formation=mother_element.actual_time_since_cells_formation,
                                                 actual_time_since_growth_stopped=mother_element.actual_time_since_growth_stopped,
                                                 actual_time_since_death=mother_element.actual_time_since_death,

                                                 thermal_time_since_primordium_formation=mother_element.thermal_time_since_primordium_formation,
                                                 thermal_time_since_emergence=mother_element.thermal_time_since_emergence,
                                                 thermal_time_since_cells_formation=mother_element.thermal_time_since_cells_formation,
                                                 thermal_potential_time_since_emergence=mother_element.thermal_potential_time_since_emergence,
                                                 thermal_time_since_growth_stopped=mother_element.thermal_time_since_growth_stopped,
                                                 thermal_time_since_death=mother_element.thermal_time_since_death,
                                                 AA=mother_element.AA,
                                                 amino_acids_consumption_by_growth_amount=0.,
                                                 amino_acids_consumption_by_growth=0.,
                                                 amino_acids_possibly_required_for_elongation=0.,
                                                 amino_acids_growth_demand=0.
                                                 )
            
            return new_child
