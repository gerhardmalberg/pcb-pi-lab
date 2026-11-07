! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Gerhard Malberg
! Elementwise coefficients, indexed by native Elmer volume element IDs.
FUNCTION PiSigma(Model, Node, Time) RESULT(Value)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Time, Value
  REAL(KIND=dp), ALLOCATABLE, SAVE :: Values(:)
  TYPE(Element_t), POINTER :: Element
  INTEGER :: N, Unit, I
  IF (.NOT. ALLOCATED(Values)) THEN
    OPEN(NEWUNIT=Unit, FILE='sigma.dat', STATUS='old', ACTION='read')
    READ(Unit,*) N
    ALLOCATE(Values(N))
    READ(Unit,*) Values
    CLOSE(Unit)
  END IF
  Element => GetCurrentElement()
  I = Element % ElementIndex
  IF (I < 1 .OR. I > SIZE(Values)) CALL Fatal('PiSigma','Invalid volume element index')
  Value = Values(I)
END FUNCTION PiSigma

FUNCTION PiHeat(Model, Node, Time) RESULT(Value)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Time, Value
  REAL(KIND=dp), ALLOCATABLE, SAVE :: Values(:)
  TYPE(Element_t), POINTER :: Element
  INTEGER :: N, Unit, I
  IF (.NOT. ALLOCATED(Values)) THEN
    OPEN(NEWUNIT=Unit, FILE='heat.dat', STATUS='old', ACTION='read')
    READ(Unit,*) N
    ALLOCATE(Values(N))
    READ(Unit,*) Values
    CLOSE(Unit)
  END IF
  Element => GetCurrentElement()
  I = Element % ElementIndex
  IF (I < 1 .OR. I > SIZE(Values)) CALL Fatal('PiHeat','Invalid volume element index')
  Value = Values(I)
END FUNCTION PiHeat
