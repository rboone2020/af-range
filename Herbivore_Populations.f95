subroutine Herbivore_Dynamics (icell)
!**** Simulate herbviore population dynamics
!****
!**** Herbivore populations (up to 10) will be represented in African landscape cells.  Each cell will hold populations.
!**** The populations will not need to track male/female, given the large numbers involved.  Nor will age-classes be tracked.
!**** Forage quality and quantity will be the primary driver of population change.  Population dynamics will come from
!**** my favored pathway of modeling dynamics (e.g., expected body mass versus actual yielding condition indices, indices
!**** affecting birth and death rates) and informed by Dangal et al. (2016), "Integrating hervibore population dynamics
!**** into a global land biosphere model: plugging animals into the Earth system" in Journal of Advances in Modeling Earth
!**** Systems. As odd as it is given my agent-based background, here all animals will be the same size, all adult, and
!**** sexes won't be tracked.
!****
!**** From STRUCTURES:
!**** integer    :: herbivore(MAX_HERBIVORES)             ! The number of herbivore species represented in the model
!**** real       :: prop_forage_avail                     ! The proportion of forage available for herbivores, which will
!****                                                     ! be cell-specific to capture the many potential differences between
!****                                                     ! food-based capacity and numbers of animals actually on the landscape (e.g., access)
!****
!**** R. Boone   Last modified: March 21, 2025    Herbivore populations were capped at a given density, a multiple of the observed density.
!                               March 5, 2026     Maximum populations made species specific.
  use Parameter_Vars
  use Structures
  implicit none
  integer  :: icell

  call Herbivore_Grazing (icell)                          ! Yields grazing fraction, reduced biomass in pools, fme for mean animal/spp/cell
  call Herbivore_Energy_Used (icell)                      ! Yields average FME acquired per animal of a given species
  call Herbivore_Weight_Change (icell)                    ! Yields a change in Herbivore_biomass(ispp) associated with the animals' condition
  call Herbivore_Mortality (icell)                        ! May be disabled.  Changes population deaths with the animals' condition
  call Herbivore_Give_Birth (icell)                       ! May be disabled.  Changes population births with the animals' condition
  call Herbivore_Pop_Limit (icell)                        ! May be disabled.  Checks if populations exceed a set limit.

end subroutine


subroutine Herbivore_Grazing (icell)
   ! Simulate grazing of animals. Distribute grazing among the pools of biomass
   ! and calculate a grazed fraction for use by L-Range and the typical Century logic.
   ! Note that in the following an if-else structure may be used to ensure that the last
   ! increment of forage is eaten, but it is for the entire 100 km2 or 400 km2 cell, and
   ! so a fractionally small amount of forage.  No need to slow the effort down with
   ! many if-elses. The forage increment is going to be larger than the ungrazed
   ! portion on an entire landscape cell.  Note: Now average per animal FME.
   !
   ! Yields FME (mean MJ/Animal for each species in cell) and grazing fraction (e.g., 0.025 offtake for month).
   !
   ! R. Boone.  Last changed:  RBB  June 20, 2023.  Repairing grazed_fraction.
   !
   use Parameter_Vars
   use Structures
   implicit none
   integer :: ispp, ipool, icell, itemp, attempt
   integer :: days_in_period
   real    :: average_kg, forage_increment
   real    :: offtake, max_pool_offtake_kgs, total_offtake_kgs, daily_offtake_percent, temp_offtake
   real    :: rand, waste, metabolizability, total_biomass, grazed_biomass, actual_offtake_kgs(8)
   real    :: sum_sought, ungrazable_green_kgs, ungrazable_dead_kgs, sum_offtake
   real    :: area_available, area_sought, area_offtake, temp_avail
   logical :: sufficient_biomass

   call Zero_Out_Grazing_Fractions (icell)
   actual_offtake_kgs = 0.0                                            ! Clearing an array
                                                                       ! Available biomass is being included in the following ...
   call Store_Temporary_Biomasses (icell)                              ! g / m^2 converted to kg / cell here
   days_in_period = MONTH_DAYS(month)
                                                                       ! If this is needlessly slow, increase the increment.  Decrease it for more resolution.
   do ispp=1,herbivores
     sufficient_biomass = .TRUE.
     call Check_for_Sufficient_Biomass (icell, ispp, sufficient_biomass)    ! Independent of the availability concept, don't let animals graze to absolutely bare ground.
     if (sufficient_biomass .eqv. .TRUE.) then                              ! There is enough on the ground not to eat it to absolute bare
       if (Rng(icell)%herbivore(ispp) .gt. 0) then                           ! If there are no animals, skip the following
         total_offtake_kgs = Rng(icell)%herbivore_biomass(ispp) * days_in_period * Herbivore_Parms%daily_offtake_prop(ispp) * &
                             Rng(icell)%herbivore(ispp)
         forage_increment = total_offtake_kgs * 0.005                        ! Going to feed animals in increments to spread consumption across pools.  This may be too slow, but could distribute forage fairly but not stringently.
         if (forage_increment .lt. 5.0) then                                 ! If there is a tiny (5 kg) amount of forage needed, set it to the increment
           forage_increment = total_offtake_kgs
         end if
         actual_offtake_kgs = 0.0                                            ! Array cleared
         attempt = 0
         do while (total_offtake_kgs .gt. forage_increment .and. attempt .lt. 250)
           call random_number(rand)
           ipool = int(rand * BIO_POOLS) + 1
           call random_number(rand)
           attempt = attempt + 1
           if (rand .le. Herbivore_Parms%pool_preference(ispp, ipool)) then                   ! Assess the preference for the biomass pool
             select case (ipool)
               case (1)
                 if (Rng(icell)%tgreen_herb_biomass .gt. forage_increment) then
                   Rng(icell)%tgreen_herb_biomass = Rng(icell)%tgreen_herb_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(1) = actual_offtake_kgs(1) + forage_increment
                   !write(*,*) 'GRASS  OFFTAKE: ', ispp, icell, Rng(icell)%tgreen_herb_biomass, actual_offtake_kgs(1)
                   attempt = 0
                 end if
               case (2)
                 if (Rng(icell)%tdead_herb_biomass .gt. forage_increment) then
                   Rng(icell)%tdead_herb_biomass = Rng(icell)%tdead_herb_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(2) = actual_offtake_kgs(2) + forage_increment
                   !write(*,*) 'DEAD GRASS  OFFTAKE: ', ispp, icell, Rng(icell)%tdead_herb_biomass, actual_offtake_kgs(2)
                   attempt = 0
                 end if
               case (3)
                 if (Rng(icell)%tgreen_shrub_biomass .gt. forage_increment) then
                   Rng(icell)%tgreen_shrub_biomass = Rng(icell)%tgreen_shrub_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(3) = actual_offtake_kgs(3) + forage_increment
                   !write(*,*) 'SHRUB  OFFTAKE: ', ispp, icell, Rng(icell)%tgreen_shrub_biomass, actual_offtake_kgs(3)
                   attempt = 0
                 end if
               case (4)
                 if (Rng(icell)%tdead_shrub_biomass .gt. forage_increment) then
                   Rng(icell)%tdead_shrub_biomass = Rng(icell)%tdead_shrub_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(4) = actual_offtake_kgs(4) + forage_increment
                   !write(*,*) 'DEAD SHRUB  OFFTAKE: ', ispp, icell, Rng(icell)%tdead_shrub_biomass, actual_offtake_kgs(4)
                   attempt = 0
                 end if
               case (5)
                 if (Rng(icell)%tfbranch_shrub_biomass .gt. forage_increment) then
                   Rng(icell)%tfbranch_shrub_biomass = Rng(icell)%tfbranch_shrub_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(5) = actual_offtake_kgs(5) + forage_increment
                   !write(*,*) 'FB SHRUB  OFFTAKE: ', ispp, icell, Rng(icell)%tfbranch_shrub_biomass, actual_offtake_kgs(5)
                   attempt = 0
                 end if
               case (6)
                 if (Rng(icell)%tgreen_tree_biomass .gt. forage_increment) then
                   Rng(icell)%tgreen_tree_biomass = Rng(icell)%tgreen_tree_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(6) = actual_offtake_kgs(6) + forage_increment
                   !write(*,*) 'TREE  OFFTAKE: ', ispp, icell, Rng(icell)%tgreen_tree_biomass, actual_offtake_kgs(6)
                   attempt = 0
                 end if
               case (7)
                 if (Rng(icell)%tdead_tree_biomass .gt. forage_increment) then
                   Rng(icell)%tdead_tree_biomass = Rng(icell)%tdead_tree_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(7) = actual_offtake_kgs(7) + forage_increment
                   !write(*,*) 'DEAD TREE  OFFTAKE: ', ispp, icell, Rng(icell)%tdead_tree_biomass, actual_offtake_kgs(7)
                   attempt = 0
                 end if
               case (8)
                 if (Rng(icell)%tfbranch_tree_biomass .gt. forage_increment) then
                   Rng(icell)%tfbranch_tree_biomass = Rng(icell)%tfbranch_tree_biomass - forage_increment
                   total_offtake_kgs = total_offtake_kgs - forage_increment
                   actual_offtake_kgs(8) = actual_offtake_kgs(8) + forage_increment
                   !write(*,*) 'FB TREE  OFFTAKE: ', ispp, icell, Rng(icell)%tfbranch_tree_biomass, actual_offtake_kgs(8)
                   attempt = 0
                 end if
             end select
           end if
         end do
         ! Convert material grazed to fermentable metabolizable energy.
         ! Gross energy content of plant tissues (MJ/kg) (18.5 from McDonald et al. 1988) As used in Savanna
         Rng(icell)%fme(ispp) = 0.0                                                   ! Array assignment. NOTE FME is the total for the entire herd in the landscape cell at this point in the code.
         do ipool=1,BIO_POOLS
           waste = actual_offtake_kgs(ipool) * Herbivore_Parms%wasted_prop(ipool)
          ! ZZZZZZ
        !if (icell .eq. 8888 .and. ispp .eq. 1) then
        ! write(*,*) 'Grazing  POPULATION : ', icell, ispp, ipool, Rng(icell)%herbivore(ispp)
        ! write(*,*) 'Grazing  FME        : ', icell, ispp, ipool, Rng(icell)%fme(ispp)
        ! write(*,*) 'Grazing  AOK        : ', icell, ispp, ipool, actual_offtake_kgs(ipool)
        ! write(*,*) 'Grazing  WASTE PROP : ', icell, ispp, ipool, Herbivore_Parms%wasted_prop(ipool)
        ! write(*,*) 'Grazing  WASTE      : ', icell, ispp, ipool, waste
        ! write(*,*) 'Grazing  DIGEST POOL: ', icell, ispp, ipool, Herbivore_Parms%digestibility_fraction(ispp, ipool)
        ! write(*,*) 'Grazing  META       : ', icell, ispp, ipool, Herbivore_Parms%metabilizability(ispp)
        ! write(*,*) 'Grazing               '
        !end if
           Rng(icell)%fme(ispp) = Rng(icell)%fme(ispp) + (( actual_offtake_kgs(ipool) - waste ) * &
                              Herbivore_Parms%digestibility_fraction(ispp, ipool) * &
                              Herbivore_Parms%metabilizability(ispp) * 18.5)
        !if (icell .eq. 8888 .and. ispp .eq. 1) then
        ! write(*,*) 'Grazing AFTER FME   : ', icell, ispp, ipool, Rng(icell)%fme(ispp)
        ! write(*,*) 'Grazing               '
        !end if
         end do

         Rng(icell)%fme(ispp) = Rng(icell)%fme(ispp) / Rng(icell)%herbivore(ispp)     ! (Divide by zero accounted for above.)
         Rng(icell)%fme(ispp) = Rng(icell)%fme(ispp) / days_in_period
        !if (icell .eq. 8888 .and. ispp .eq. 1) then
        ! write(*,*) 'Grazing ENERGY ACQUIRED: ', icell, ispp, Rng(icell)%fme(ispp)
        ! write(*,*) 'Grazing                 '
        !end if
       end if                                                                         ! No need to set fractions to 0 if no animals ... already done by zero_out_grazing_fraction above
     end if
   end do
   call Get_Grazed_Fraction (icell, actual_offtake_kgs)

end subroutine


subroutine Zero_Out_Grazing_Fractions (icell)
  ! Zero-out the grazing fraction for the cell.
  ! R.Boone  Last edit:  April 17, 2023
  use Parameter_Vars
  use Structures
  implicit none
  integer :: icell

  Rng(icell)%grazed_green_herb_fraction = 0.0
  Rng(icell)%grazed_dead_herb_fraction = 0.0
  Rng(icell)%grazed_green_shrub_fraction = 0.0
  Rng(icell)%grazed_dead_shrub_fraction = 0.0
  Rng(icell)%grazed_fbranch_shrub_fraction = 0.0
  Rng(icell)%grazed_green_tree_fraction = 0.0
  Rng(icell)%grazed_dead_tree_fraction = 0.0
  Rng(icell)%grazed_fbranch_tree_fraction = 0.0

end subroutine


subroutine Store_Temporary_Biomasses (icell)
  ! Store the biomasses of cells in temporary variables that will be decremented to 0 as grazing is simulated.
  ! *** Note that the proportion of vegetation available as forage is incorporated here ***
  ! R.Boone  Last edit:  February 21, 2025.  Altering forage calculations to improve understory estimates.
  use Parameter_Vars
  use Structures
  implicit none

  integer :: icell, iunit
  real    :: available, cell_res
  real    :: h_frac, s_frac, t_frac
  real    :: h_green_biomass, h_dead_biomass
  real    :: s_green_biomass, s_dead_biomass, s_fb_biomass
  real    :: t_green_biomass, t_dead_biomass, t_fb_biomass

  available = (Area(Rng(icell)%x, Rng(icell)%y)%prop_forage_avail / 100.) / 12.0      ! NOTE: Conversion to monthly value.  Could have month-specific.  Converts from % to proportion.
  cell_res = Herbivore_Parms%cell_resolution

  iunit = Rng(icell)%range_type

  ! LEAF_CARBON(*), for example, stores the carbon for that plant type, regardless of whether it is in the overstory or understory.  If plants covered the entire area,
  ! it would be simply LEAF_CARBON plus other components, but they don't have full cover.  Most straightforward (and using the available data) would be to compare
  ! populations to the total potential population for the reference 1 km2 area.   (potential populations can't be 0, so division should not need trapped)
  h_frac = (Rng(icell)%total_population(H_LYR) + Rng(icell)%total_population(H_S_LYR) + &
            Rng(icell)%total_population(H_T_LYR)) / Parms(iunit)%pot_population(H_FACET)
  s_frac = (Rng(icell)%total_population(S_LYR) + Rng(icell)%total_population(S_T_LYR)) / Parms(iunit)%pot_population(S_FACET)
  t_frac =  Rng(icell)%total_population(T_LYR) / Parms(iunit)%pot_population(H_FACET)                                                      ! Equal to facet

  ! BIOMASS COMPONENTS  Note: Herbs, shrubs, and trees are summed eventually
  h_green_biomass = ((Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET)) * 2.5) * h_frac
  h_dead_biomass  = ((Rng(icell)%dead_leaf_carbon(H_FACET) + Rng(icell)%dead_seed_carbon(H_FACET)) * 2.5) * h_frac
  s_green_biomass = ((Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET)) * 2.5) * s_frac
  s_dead_biomass =  ((Rng(icell)%dead_leaf_carbon(S_FACET) + Rng(icell)%dead_seed_carbon(S_FACET)) * 2.5) * s_frac
  s_fb_biomass   =  (Rng(icell)%fine_branch_carbon(S_FACET) * 2.5) * s_frac
  t_green_biomass = ((Rng(icell)%leaf_carbon(T_FACET) + Rng(icell)%seed_carbon(T_FACET)) * 2.5) * t_frac
  t_dead_biomass =  ((Rng(icell)%dead_leaf_carbon(T_FACET) + Rng(icell)%dead_seed_carbon(T_FACET)) * 2.5) * t_frac
  t_fb_biomass =    (Rng(icell)%fine_branch_carbon(T_FACET) * 2.5) * t_frac

  ! NOTE Temporary storage locations are in KG per cell
  Rng(icell)%tgreen_herb_biomass  = (h_green_biomass * available * cell_res * cell_res) / 1000.0
  Rng(icell)%tdead_herb_biomass   = (h_dead_biomass  * available * cell_res * cell_res) / 1000.0
  Rng(icell)%tgreen_shrub_biomass = (s_green_biomass * available * cell_res * cell_res) / 1000.0
  Rng(icell)%tdead_shrub_biomass  = (s_dead_biomass * available * cell_res * cell_res)  / 1000.0
  Rng(icell)%tfbranch_shrub_biomass = (s_fb_biomass * available * cell_res * cell_res)  / 1000.0
  Rng(icell)%tgreen_tree_biomass = (t_green_biomass * available * cell_res * cell_res) / 1000.0
  Rng(icell)%tdead_tree_biomass  = (t_dead_biomass * available * cell_res * cell_res)  / 1000.0
  Rng(icell)%tfbranch_tree_biomass = (t_fb_biomass * available * cell_res * cell_res)  / 1000.0

 !if (icell .eq. 8888) then
 !  write(*,*) 'BIOMASSES  icell, gherb, dherb, gshrub, dshrub, fbshrub, gtree, dtree, fbtree '
 !  write(*,*) icell, Rng(icell)%tgreen_herb_biomass, Rng(icell)%tdead_herb_biomass, Rng(icell)%tgreen_shrub_biomass, &
 !   Rng(icell)%tdead_shrub_biomass, Rng(icell)%tfbranch_shrub_biomass, Rng(icell)%tgreen_tree_biomass, &
 !   Rng(icell)%tdead_tree_biomass, Rng(icell)%tfbranch_tree_biomass
 !  write(*,*) 'HERBIVORES icell, cattle, goats, sheep, camels. donkeys, grazers, mixed '
 !  write(*,*) icell, Rng(icell)%herbivore(1), Rng(icell)%herbivore(2), Rng(icell)%herbivore(3), Rng(icell)%herbivore(4), &
 !    Rng(icell)%herbivore(5), Rng(icell)%herbivore(6), Rng(icell)%herbivore(7)
 !end if

end subroutine


subroutine Check_for_Sufficient_Biomass (icell, ispp, sufficient_biomass)
  ! Regardless of the idea that some forage will be unavailable, we must prevent the animals from eating the
  ! forage all the way to bare ground.  They won't do that in general.
  ! No need to convert to cell biomasses here, in so far as those would be scaler changes that wouldn't affect logic.
  ! This will favor green or dead shrub and trees.  Leave them out and focus on herbs only?  Or include them all?
  ! Here I will include them all. A cow that has no grass to eat may nibble on green tree biomass.
  ! Putting branches with dead material.
  !
  ! R.Boone  Last edit:  February 21, 2025
  use Parameter_Vars
  use Structures
  implicit none

  integer :: icell, ispp, iunit
  real    :: green_biomass, dead_biomass
  real    :: h_green_biomass, h_dead_biomass, s_green_biomass, s_dead_biomass, t_green_biomass, t_dead_biomass
  logical :: sufficient_biomass
  real    :: h_frac, s_frac, t_frac

  sufficient_biomass = .TRUE.

  green_biomass = 0.0
  dead_biomass = 0.0

  iunit = Rng(icell)%range_type

  ! LEAF_CARBON(*), for example, stores the carbon for that plant type, regardless of whether it is in the overstory or understory.  If plants covered the entire area,
  ! it would be simply LEAF_CARBON plus other components, but they don't have full cover.  Most straightforward (and using the available data) would be to compare
  ! populations to the total potential population for the reference 1 km2 area.   (potential populations can't be 0, so division should not need trapped)
  h_frac = (Rng(icell)%total_population(H_LYR) + Rng(icell)%total_population(H_S_LYR) + &
            Rng(icell)%total_population(H_T_LYR)) / Parms(iunit)%pot_population(H_FACET)
  s_frac = (Rng(icell)%total_population(S_LYR) + Rng(icell)%total_population(S_T_LYR)) / Parms(iunit)%pot_population(S_FACET)
  t_frac =  Rng(icell)%total_population(T_LYR) / Parms(iunit)%pot_population(H_FACET)                                                      ! Equal to facet

  ! BIOMASS COMPONENTS  Note: Herbs, shrubs, and trees are summed eventually
  h_green_biomass = ((Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET)) * 2.5) * h_frac
  h_dead_biomass  = ((Rng(icell)%dead_leaf_carbon(H_FACET) + Rng(icell)%dead_seed_carbon(H_FACET)) * 2.5) * h_frac
  s_green_biomass = ((Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET)) * 2.5) * s_frac
  s_dead_biomass =  ((Rng(icell)%dead_leaf_carbon(S_FACET) + Rng(icell)%dead_seed_carbon(S_FACET)) * 2.5) * s_frac
  s_dead_biomass =  s_dead_biomass + (Rng(icell)%fine_branch_carbon(S_FACET) * 2.5) * s_frac
  t_green_biomass = ((Rng(icell)%leaf_carbon(T_FACET) + Rng(icell)%seed_carbon(T_FACET)) * 2.5) * t_frac
  t_dead_biomass =  ((Rng(icell)%dead_leaf_carbon(T_FACET) + Rng(icell)%dead_seed_carbon(T_FACET)) * 2.5) * t_frac
  t_dead_biomass =  t_dead_biomass + (Rng(icell)%fine_branch_carbon(T_FACET) * 2.5) * t_frac

  green_biomass = h_green_biomass + s_green_biomass + t_green_biomass
  dead_biomass  = h_dead_biomass  + s_dead_biomass  + t_dead_biomass

  ! Stop the foraging if green is gone.  If green is gone but some dead remains, don't stop foraging.
  if (green_biomass .lt. Herbivore_Parms%ungrazable_biomass(ispp, 1)) then
    sufficient_biomass = .FALSE.
    if (dead_biomass .ge. Herbivore_Parms%ungrazable_biomass(ispp, 2)) then
      sufficient_biomass = .TRUE.
    end if
  end if
  ! YYYY
 !if (icell .eq. 8888 .and. ispp .eq. 1) then
 !  write(*,*) 'LEAF SEED H_FRAC   ', icell, Rng(icell)%leaf_carbon(H_FACET), Rng(icell)%seed_carbon(H_FACET), &
 !   h_frac, Rng(icell)%facet_cover(H_FACET)
 !  write(*,*) 'LEAF SEED S_FRAC   ', icell, Rng(icell)%leaf_carbon(S_FACET), Rng(icell)%seed_carbon(S_FACET), s_frac, &
 !   Rng(icell)%facet_cover(S_FACET)
 !  write(*,*) 'LEAF SEED T_FRAC   ', icell, Rng(icell)%leaf_carbon(T_FACET), Rng(icell)%seed_carbon(T_FACET), t_frac, &
 !   Rng(icell)%facet_cover(T_FACET)
 !
 !  write(*,*) 'CHECKING  icell, green, dead, sufficient ', icell, green_biomass, dead_biomass, sufficient_biomass
 !end if
end subroutine


subroutine Get_Grazed_Fraction (icell, offtake_kgs)
  ! Store the fraction of total biomass consumed in the cell.  This involves recalculation, but I don't want to store yet
  ! eight more layers for all of Africa. NOTE: Yields proportion for each biomass pool.
  ! R.Boone  Last edit:  June 20, 2023   Offtake handled more explicitly using passed array.
  use Parameter_Vars
  use Structures
  implicit none

  integer :: icell, iunit
  real    :: cell_res, temp_biomass, offtake_kgs(8)
  real    :: h_green_biomass, h_dead_biomass
  real    :: s_green_biomass, s_dead_biomass, s_fb_biomass
  real    :: t_green_biomass, t_dead_biomass, t_fb_biomass
  real    :: h_frac, s_frac, t_frac
  real    :: tgreen_herb_bio, tdead_herb_bio
  real    :: tgreen_shrub_bio, tdead_shrub_bio, tfbranch_shrub_bio
  real    :: tgreen_tree_bio,  tdead_tree_bio,  tfbranch_tree_bio

  cell_res = Herbivore_Parms%cell_resolution

  iunit = Rng(icell)%range_type

  h_frac = (Rng(icell)%total_population(H_LYR) + Rng(icell)%total_population(H_S_LYR) + &
            Rng(icell)%total_population(H_T_LYR)) / Parms(iunit)%pot_population(H_FACET)
  s_frac = (Rng(icell)%total_population(S_LYR) + Rng(icell)%total_population(S_T_LYR)) / Parms(iunit)%pot_population(S_FACET)
  t_frac =  Rng(icell)%total_population(T_LYR) / Parms(iunit)%pot_population(H_FACET)                                                      ! Equal to facet

  h_green_biomass = ((Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET)) * 2.5) * h_frac
  h_dead_biomass  = ((Rng(icell)%dead_leaf_carbon(H_FACET) + Rng(icell)%dead_seed_carbon(H_FACET)) * 2.5) * h_frac
  s_green_biomass = ((Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET)) * 2.5) * s_frac
  s_dead_biomass  = ((Rng(icell)%dead_leaf_carbon(S_FACET) + Rng(icell)%dead_seed_carbon(S_FACET)) * 2.5) * s_frac
  s_fb_biomass    = (Rng(icell)%fine_branch_carbon(S_FACET) * 2.5) * s_frac
  t_green_biomass = ((Rng(icell)%leaf_carbon(T_FACET) + Rng(icell)%seed_carbon(T_FACET)) * 2.5) * t_frac
  t_dead_biomass  = ((Rng(icell)%dead_leaf_carbon(T_FACET) + Rng(icell)%dead_seed_carbon(T_FACET)) * 2.5) * t_frac
  t_fb_biomass    = (Rng(icell)%fine_branch_carbon(T_FACET) * 2.5) * t_frac

  tgreen_herb_bio    = (h_green_biomass * cell_res * cell_res) / 1000.0                     ! "available" not used here, it is relative to total production
  tdead_herb_bio     = (h_dead_biomass  * cell_res * cell_res) / 1000.0
  tgreen_shrub_bio   = (s_green_biomass * cell_res * cell_res) / 1000.0
  tdead_shrub_bio    = (s_dead_biomass * cell_res * cell_res)  / 1000.0
  tfbranch_shrub_bio = (s_fb_biomass * cell_res * cell_res)  / 1000.0
  tgreen_tree_bio    = (t_green_biomass * cell_res * cell_res) / 1000.0
  tdead_tree_bio     = (t_dead_biomass * cell_res * cell_res)  / 1000.0
  tfbranch_tree_bio  = (t_fb_biomass * cell_res * cell_res)  / 1000.0

  Rng(icell)%grazed_green_herb_fraction = 0.0
  Rng(icell)%grazed_dead_herb_fraction = 0.0
  Rng(icell)%grazed_green_shrub_fraction = 0.0
  Rng(icell)%grazed_dead_shrub_fraction = 0.0
  Rng(icell)%grazed_fbranch_shrub_fraction = 0.0
  Rng(icell)%grazed_green_tree_fraction = 0.0
  Rng(icell)%grazed_dead_tree_fraction = 0.0
  Rng(icell)%grazed_fbranch_tree_fraction = 0.0

  if (offtake_kgs(1) .gt. 0.01 .and. tgreen_herb_bio .gt. 0.01) then
    Rng(icell)%grazed_green_herb_fraction = offtake_kgs(1) / tgreen_herb_bio
  end if
  if (offtake_kgs(2) .gt. 0.01 .and. tdead_herb_bio .gt. 0.01) then
    Rng(icell)%grazed_dead_herb_fraction = offtake_kgs(2) / tdead_herb_bio
  end if
  if (offtake_kgs(3) .gt. 0.01 .and. tgreen_shrub_bio .gt. 0.01) then
    Rng(icell)%grazed_green_shrub_fraction = offtake_kgs(3) / tgreen_shrub_bio
  end if
  if (offtake_kgs(4) .gt. 0.01 .and. tdead_shrub_bio .gt. 0.01) then
    Rng(icell)%grazed_dead_shrub_fraction = offtake_kgs(4) / tgreen_shrub_bio
  end if
  if (offtake_kgs(5) .gt. 0.01 .and. tfbranch_shrub_bio .gt. 0.01) then
    Rng(icell)%grazed_fbranch_shrub_fraction = offtake_kgs(5) / tfbranch_shrub_bio
  end if
  if (offtake_kgs(6) .gt. 0.01 .and. tgreen_tree_bio .gt. 0.01) then
    Rng(icell)%grazed_green_tree_fraction = offtake_kgs(6) / tgreen_tree_bio
  end if
  if (offtake_kgs(7) .gt. 0.01 .and. tdead_tree_bio .gt. 0.01) then
    Rng(icell)%grazed_dead_tree_fraction = offtake_kgs(7) / tdead_tree_bio
  end if
  if (offtake_kgs(8) .gt. 0.01 .and. tfbranch_tree_bio .gt. 0.01) then
    Rng(icell)%grazed_fbranch_tree_fraction = offtake_kgs(8) / tfbranch_tree_bio
  end if

  ! XXX
  !if (icell .eq. 8888) then
  !  write(*,*) 'FRACTIONS H S T      ', icell, h_frac, s_frac, t_frac
  !  write(*,*) 'Biomasses            ', icell, tgreen_herb_bio, tdead_herb_bio
  !  write(*,*) 'Offtake 1 and biomass', icell, offtake_kgs(1)
  !  write(*,*) 'FRACTION GREEN HERB  ', icell, Rng(icell)%grazed_green_herb_fraction
  !  write(*,*) 'FRACTION DEAD  HERB  ', icell, Rng(icell)%grazed_dead_herb_fraction
  !  write(*,*) 'FRACTION GREEN SHRUB ', icell, Rng(icell)%grazed_green_shrub_fraction
  !  write(*,*) 'FRACTION DEAD  SHRUB ', icell, Rng(icell)%grazed_dead_shrub_fraction
  !  write(*,*) 'FRACTION FB    SHRUB ', icell, Rng(icell)%grazed_fbranch_shrub_fraction
  !  write(*,*) 'FRACTION GREEN TREE  ', icell, Rng(icell)%grazed_green_tree_fraction
  !  write(*,*) 'FRACTION DEAD  TREE  ', icell, Rng(icell)%grazed_dead_tree_fraction
  !  write(*,*) 'FRACTION FB    TREE  ', icell, Rng(icell)%grazed_fbranch_tree_fraction
  !end if

end subroutine


subroutine Herbivore_Weight_Change (icell)
   ! Caculate weight change for animals based on their energy use and energy acquired
   !
   ! R.Boone    Last changed:  May 23, 2023
   use Parameter_Vars
   use Structures
   implicit none
   integer :: ispp, icell, ipool, days_in_period
   real    :: delta_kg  ! ratio_ne_to_me not needed here, energies are calculated as net energies already
   real    :: max_mass, min_mass

   days_in_period = MONTH_DAYS(month)

   do ispp=1,herbivores
     if (Rng(icell)%herbivore(ispp) .eq. 0) then                             ! No animals of that species, so no average biomass
       Rng(icell)%herbivore_biomass(ispp) = 0
     else
       if (Rng(icell)%fme(ispp) .gt. (Rng(icell)%energy_used(ispp) * days_in_period)) then
         ! Animals will gain weight.  Live weight increased ... 26 MJ NE/kg live weight
         ! 26 MJ NE/kg live weight is the contribution required to gain a kg.                                                                             (Coppock et al. 1983, ARC 1980)
         ! Weight gain must be limited here, as an average response, as well as below.
         ! FME_ACQUIRED and ENERGY_USED are in MJ/d/An.  Below the DELTA_KG is changed to a monthly weight change.
         delta_kg = ((Rng(icell)%fme(ispp) - Rng(icell)%energy_used(ispp)) * days_in_period) / 26.                  ! Daily weight gain
         delta_kg = delta_kg / days_in_period
         if (delta_kg .gt. Herbivore_Parms%min_max_mass_change(ispp, 2)) then
           delta_kg = Herbivore_Parms%min_max_mass_change(ispp, 2)
         end if
       else
         ! Animals will loose weight (or stay the same if exactly equal)
         ! Effect of body mass conversion to NE is 0.84, from Coppock.
         ! Weight loss must be limited here, as an average response, as well as below.
         ! FME_ACQUIRED and ENERGY_USED are in MJ/d/An.  Below the DELTA_KG is changed to a monthly weight change.
         delta_kg = ((Rng(icell)%fme(ispp) - Rng(icell)%energy_used(ispp)) * days_in_period) / (26.0 * 0.84)                                           ! Daily weight lost
         delta_kg = delta_kg / days_in_period
         if (delta_kg .lt. Herbivore_Parms%min_max_mass_change(ispp, 1)) then
           delta_kg = Herbivore_Parms%min_max_mass_change(ispp, 1)
         end if
       end if
       delta_kg = delta_kg * (365./12.)                                       ! Convert the daily figure to a monthly figure.
       Rng(icell)%herbivore_biomass(ispp) = Rng(icell)%herbivore_biomass(ispp) + delta_kg

       ! Delta_kg may be positive or negative.
       ! DO NOT ALLOW THE COHORT OF ANIMALS TO GET TOO HEAVY (> 1.0 CONDITION INDEX)
       ! OR TOO LIGHT (< 0.0 CONDITION INDEX).  Can't allow mouse size or elephant size cattle, for example.
       min_mass = Herbivore_Parms%expected_kg(ispp) * Herbivore_Parms%min_max_mass_ratio(ispp, 1)
       max_mass = Herbivore_Parms%expected_kg(ispp) * Herbivore_Parms%min_max_mass_ratio(ispp, 2)
       if (Rng(icell)%herbivore_biomass(ispp) .gt. max_mass) then
         Rng(icell)%herbivore_biomass(ispp) = max_mass
       end if
       if (Rng(icell)%herbivore_biomass(ispp) .lt. min_mass) then
         Rng(icell)%herbivore_biomass(ispp) = min_mass
       end if                                                         ! end if, fme > energy used

     !  ! WWW
     !  if (icell .eq. 8888 .and. ispp .eq. 1) then
     !    write(*,*) 'WEIGHT CHANGE  icell, ispp, acquired, used, mass, pop: ', icell, ispp, &
     !    Rng(icell)%fme(ispp), Rng(icell)%energy_used(ispp), Rng(icell)%herbivore_biomass(ispp), Rng(icell)%herbivore(ispp)
     !  end if
     end if                                                           ! end if, 0 animals or not
   end do

end subroutine


subroutine Herbivore_Energy_Used (icell)
   ! Caculate how much energy animals are using.  Little can be included for energy use on large cells
   ! and single populations for each herd.  An exception is thermal stress, which is regionally applicable.
   ! Voluntary energy may be represented as well.
   !
   ! Yields: The outcome from this procedure is per day, mean MJ per animal.   NOTE: Per day.
   !
   ! R. Boone   Last modified:  May 16, 2023
   use Parameter_Vars
   use Structures
   implicit none

   integer  :: ispp, i, icell
   real     :: average_kg, average_ci, effect_condition, effect_thermal, base_energy_used
   real     :: vol_energy, alint, rand, expected_kg, b_range, voluntary_energy_use(6)
   real     :: average_temperature, condition_fraction

   do ispp=1, herbivores
     if (Rng(icell)%herbivore(ispp) .gt. 0) then
       average_kg = Rng(icell)%herbivore_biomass(ispp)
       expected_kg = Herbivore_Parms%expected_kg(ispp)
       ! Baseline energy use is the minimum basal metabolism times the average body mass of animals in the cell.
       base_energy_used = Herbivore_Parms%min_max_energy_use(ispp, 1) * average_kg                                          ! MJ/d/Animal
       b_range = Herbivore_Parms%min_max_mass_ratio(ispp, 2) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)
       average_ci = ((average_kg / expected_kg) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)) / b_range

       ! Effect of condition on voluntary energy use.  Animals that are in poor condition can cut-back on energy needs. A multiplier on basal energy use.
       do i=1,6
         voluntary_energy_use(i) = Herbivore_Parms%max_energy_use_vs_ci(ispp, i)
       end do
       condition_fraction = alint(average_ci, voluntary_energy_use, 3, 'vol_eff     ')
       effect_condition = (condition_fraction * (Herbivore_Parms%min_max_energy_use(ispp, 2) - &
                          Herbivore_Parms%min_max_energy_use(ispp, 1))) * average_kg
       ! Thermal stress.  Note that in other applications, cold is represented, but extreme heat is not.  MODIFY THIS?    Not going to use bedding thermal limit at this time, since the simpler form of energy use is represented.
       average_temperature = (Area(Rng(icell)%x,Rng(icell)%y)%max_temp + Area(Rng(icell)%x,Rng(icell)%y)%min_temp) / 2.0
       if (average_temperature .le. Herbivore_Parms%critical_temps(ispp, 1)) then
         effect_thermal = Herbivore_Parms%thermal_cost_prop(ispp) * (average_temperature - Herbivore_Parms%critical_temps(ispp, 1))
       else
         effect_thermal = 0.                                                                                                 ! MJ/d/Animal
       end if
       effect_thermal = effect_thermal * average_kg
       Rng(icell)%energy_used(ispp) = base_energy_used + effect_thermal + effect_condition
       ! Make sure energy used doesn't exceed the maximum for the species
       if (Rng(icell)%energy_used(ispp) .gt. (Herbivore_Parms%min_max_energy_use(ispp, 2) * average_kg)) then
         Rng(icell)%energy_used(ispp) = Herbivore_Parms%min_max_energy_use(ispp, 2) * average_kg                               ! MJ/d/Animal
       end if
     else
       Rng(icell)%energy_used(ispp) = 0.0
     end if

   !  ! VVVV
   !  if (icell .eq. 8888 .and. ispp .eq. 1) then
   !    write(*,*) 'ENERGY USED:  icell, ispp, energy: ', icell, ispp, Rng(icell)%energy_used(ispp)
   !  end if
  end do

end subroutine


subroutine Herbivore_Mortality (icell)
   ! Simulate herbivore population decrease.  Individual dynamics are not being represented.  Individuals are not tracked
   ! and do not live and die.  Here, entire groups of animals will be removed in-line with the condition index of
   ! the population.
   !
   ! R. Boone   Last modified:  May 17, 2023
   use Parameter_Vars
   use Structures
   implicit none

   integer  :: ispp, i, icell
   real     :: average_kg, average_ci, effect_condition
   real     :: alint, rand, expected_kg, b_range, death_v_ci(8)

   do ispp=1, herbivores
     if (Herbivore_Parms%model_dynamics(ispp) .eq. 2) then
       !write(*,*) 'Herbivore Mortality'
       if (Rng(icell)%herbivore(ispp) .gt. 0) then
         average_kg = Rng(icell)%herbivore_biomass(ispp)
         expected_kg = Herbivore_Parms%expected_kg(ispp)
         b_range = Herbivore_Parms%min_max_mass_ratio(ispp, 2) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)              ! Note, not masses but ratios
         average_ci = ((average_kg / expected_kg) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)) / b_range
         do i=1,8
           death_v_ci(i) = Herbivore_Parms%death_rate_vs_ci(ispp, i)
         end do
         effect_condition = alint(average_ci, death_v_ci, 4, 'ci_die_eff  ')
        !if (icell .eq. 8888 .and. ispp .eq. 1) then
        !  write(*,*) 'MORTALITY BEFORE:  ', icell, ispp, Rng(icell)%herbivore(ispp)
        !  write(*,*) '     Average kg :  ', icell, ispp, average_kg
        !  write(*,*) '     Expected kg:  ', icell, ispp, expected_kg
        !  write(*,*) '     B range    :  ', icell, ispp, b_range
        !  write(*,*) '     Average CI :  ', icell, ispp, average_ci
        !  write(*,*) '     Effect cond:  ', icell, ispp, effect_condition
        !  write(*,*) '     DEATHS     :  ', icell, ispp, (Rng(icell)%herbivore(ispp) * effect_condition)
        !  write(*,*) 'MORTALITY AFTER :  ', icell, ispp, Rng(icell)%herbivore(ispp) - &
        !   (Rng(icell)%herbivore(ispp) * effect_condition)
        !end if
         Rng(icell)%herbivore(ispp) = Rng(icell)%herbivore(ispp) - (Rng(icell)%herbivore(ispp) * effect_condition)
       end if
     end if
  end do

end subroutine


subroutine Herbivore_Give_Birth (icell)
   ! Simulate herbivore population increase.  Individual dynamics are not being represented.  Individuals are not tracked
   ! and do not live and die.  Here, entire groups of animals will be added in-line with the condition index of
   ! the population.
   !
   ! R. Boone   Last modified:  May 17, 2023
   use Parameter_Vars
   use Structures
   implicit none

   integer  :: ispp, i, icell
   real     :: average_kg, average_ci, effect_condition
   real     :: alint, rand, expected_kg, b_range, birth_v_ci(8)

   do ispp=1, herbivores
     if (Herbivore_Parms%model_dynamics(ispp) .eq. 2) then
       !write(*,*) 'Herbivore Give Birth'
       if (Rng(icell)%herbivore(ispp) .gt. 0) then
         average_kg = Rng(icell)%herbivore_biomass(ispp)
         expected_kg = Herbivore_Parms%expected_kg(ispp)
         b_range = Herbivore_Parms%min_max_mass_ratio(ispp, 2) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)              ! Note, not masses but ratios
         average_ci = ((average_kg / expected_kg) - Herbivore_Parms%min_max_mass_ratio(ispp, 1)) / b_range

         do i=1,8
           birth_v_ci(i) = Herbivore_Parms%birth_rate_vs_ci(ispp, i)
         end do
         effect_condition = alint(average_ci, birth_v_ci, 4, 'ci_birth_eff')
         Rng(icell)%herbivore(ispp) = Rng(icell)%herbivore(ispp) + (Rng(icell)%herbivore(ispp) * effect_condition)
       ! if (icell .eq. 8888 .and. ispp .eq. 1) then
       !   write(*,*) 'BIRTHS  : ', icell, ispp, (Rng(icell)%herbivore(ispp) * effect_condition)
       ! end if
       end if
     end if
  end do

end subroutine


subroutine Herbivore_Pop_Limit (icell)
   ! If a cap is enabled, do not allow populations to exceed a multipler of the observed density.
   ! Among responses, this will prevent landscape units without animals from supporting many.
   !
   ! R. Boone   Last modified:  March 5, 2026.  The cap is now species specific.
   use Parameter_Vars
   use Structures
   implicit none

   integer  :: ispp, i, icell, iunit, max_pop_dens
   real     :: cell_res_km

   cell_res_km = Herbivore_Parms%cell_resolution / 1000.0

   iunit = Rng(icell)%range_type

   do ispp=1, herbivores
     if (Herbivore_Parms%model_dynamics(ispp) .eq. 2) then
       !write(*,*) 'Herbivore capping populations'
       max_pop_dens = (Herbivore_Parms%max_density_multiplier(ispp) * Parms(iunit)%herbivore_density(ispp)) * &
                     (cell_res_km * cell_res_km)
       if (Rng(icell)%herbivore(ispp) .gt. max_pop_dens) then
      ! if (icell .eq. 8888 .and. ispp .eq. 1) then
      !write(*,*) 'Limit:', icell, ispp, Herbivore_Parms%max_density_multiplier, Parms(iunit)%herbivore_density(ispp), cell_res_km
      !  write(*,*) '      ', max_pop_dens
      !end if
         Rng(icell)%herbivore(ispp) = max_pop_dens
       end if
     end if
   end do

end subroutine


function alint (x, data1d, imx, call_id)
   ! Does a linear interpolation between values
   ! ALINT is adapted from SAVANNA by Michael Coughenour
   ! It calculates a y value for an x value, given a series of points defining a (often broken) line.
   ! X is the x value, data1d is a one-dimensional array with paired entries, imx is the number of pairs.
   implicit none
   integer :: i=1, j=1, m=1, k=1, imx                       !Initialization of m failed.  Not sure why.
   real    :: data2d(2,4), data1d(imx*2)
   real    :: alint, x
   character(12) :: call_id                          ! Needed some way to recognize what was calling a failing ALINT function. General, that is, what was passing ALINT a -NaN value

!   if (imx > 4) then
!     write(*,*) 'ALINT, X, IMX, call_id :',x, imx, call_id
!   end if
   m=1
   do k=1,imx                                               ! Fill a temporary 2-D array for the next step
      do j=1,2
         data2d(j,k) = data1d(m)
         m=m+1
      end do
   end do
   if (x <= data2d(1,1)) then                               ! Don't let the Y values go below the min Y or above the max Y
      alint = data2d(2,1)
      return
   end if
   if (x >= data2d(1,imx)) then
      alint = data2d(2,imx)
      return
   end if
   do i=1,imx-1
      if (x <= data2d(1,i+1)) then
         k = i
         exit
      end if
   end do
   if (k > 3) then
     write(*,*) 'ALINT, X, IMX, call_id :',x, imx, call_id
   end if
   alint = data2d(2,k)+(data2d(2,k+1)-data2d(2,k)) / (data2d(1,k+1)-data2d(1,k))*(x-data2d(1,k))

end function
