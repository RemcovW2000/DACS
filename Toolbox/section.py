import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from Toolbox.helperclasses import boom, segment

class section:
    def __init__(self, members, Sx, Sy, leadingedge, trailingedge, A):
        self.members = members
        self.Sx, self.Sy = Sx, Sy
        self.leadingedge = leadingedge
        self.trailingedge = trailingedge
        self.A = A

        self.qs0 = None
        self.Corrected = False

        # Functions to run on initialisation:
        self.CalculateShearMoment()
        self.CalculateDeltas()
        self.Calculate_int_qb_ds_over_t()
        self.CalculateG()
        self.CalculatePrefactor()

    def CalculateShearMoment(self):
        print('----------------------------------------------')
        moment = 0
        for i, member in enumerate(self.members):
            membermoment = self.ShearMomentSegments(member.segments)
            print('moment in member ', i, ' : ', membermoment)
            moment += membermoment
        self.ShearMoment = moment
        print('----------------------------------------------')
        print('moment in section: ', moment)
        print('----------------------------------------------')
        return self.ShearMoment

    def ShearMomentSegments(self, segmentlist):
        # shear moment must be zero around any point (we'll take neutral point)
        moment = 0
        for segment in segmentlist:
            p1 = np.array(segment.p1)
            p2 = np.array(segment.p2)
            center = p1 + (p2 - p1)/2
            rx, ry = center[0], center[1]

            # force vector:
            directionvector = p2-p1
            # force by segments may not be consistent?!
            F = segment.qs * directionvector
            # print('force vector: ', F)
            Fx, Fy = F[0], F[1]
            M = ry * Fx - rx * Fy
            moment += M

        return moment

    def CalculateDeltas(self):
        if not self.leadingedge:
            frontmember = self.members[-1]

            difference = np.array(frontmember.endcoord) - np.array(frontmember.startcoord)

            arclength = np.linalg.norm(difference)
            self.deltafront = arclength / frontmember.panel.h
        if not self.trailingedge:
            backmember = self.members[1]
            difference = np.array(backmember.endcoord) - np.array(backmember.startcoord)

            arclength = np.linalg.norm(difference)
            self.deltaback = arclength / backmember.panel.h

        deltawhole = 0
        for member in self.members:
            # find arc length of member:
            segments = member.segments
            memberarclength = 0
            for segment in segments:
                # find the length of each segment and add to member arc length
                p1 = np.array(segment.p1)
                p2 = np.array(segment.p2)
                segmentlength = abs(np.linalg.norm(p2-p1))
                memberarclength += segmentlength

            member.arclength = memberarclength

            # add the delta of this member!
            deltawhole += memberarclength/member.panel.h

        self.deltawhole = deltawhole
        return

    def Calculate_int_qb_ds_over_t(self):
        productlist = []
        for member in self.members:
            productlistcurrent = [(segment.qs * np.linalg.norm(np.array(segment.p2)-np.array(segment.p1))) /member.panel.h for segment in member.segments]
            productlist += productlistcurrent

        integral = math.fsum(productlist)
        self.int_qb_ds_over_t = integral
        return self.int_qb_ds_over_t

    def CalculateG(self):
        # find the average g modulus of the cross section somehow?
        # TODO: find correct way to calculate G modulus
        # average wrt the arc length? -> weird, it should have something to do with it's contribution to the torsional stifness.
        Glist = [member.panel.Gxy * member.arclength for member in self.members]
        Arclengthlist = [member.arclength for member in self.members]
        self.G = sum(Glist)/sum(Arclengthlist)
        return self.G

    def CalculatePrefactor(self):
        self.prefactor = 1/(2 * self.G * self.A)
        return 1/(2 * self.G * self.A)

    def ShearCorrection(self):
        if not self.Corrected:
            for member in self.members:
                for segment in member.segments:
                    segment.qs += self.qs0
            self.Corrected = True
        else:
            print('WARNING, cross section has already been corrected! Correction not carried out.')
        return

    def PlotShearFlow(self):
        allsegments = []
        for member in self.members:
            newsegments = [(segment.p1, segment.p2, segment.qs) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        # Extract values to normalize
        values = [value for _, _, value in allsegments]
        min_value = min(values)
        max_value = max(values)
        max_abs_value = max(abs(min_value), abs(max_value))

        # Determine the appropriate colormap and normalization
        if max_value <= 0:
            # All values are negative or zero, use absolute values for intensity
            cmap = plt.get_cmap('Blues')
            norm = Normalize(vmin=0, vmax=abs(min_value))
            normalized_values = [abs(value) for value in values]
        elif min_value >= 0:
            # All values are positive or zero
            cmap = plt.get_cmap('Reds')
            norm = Normalize(vmin=0, vmax=max_value)
            normalized_values = values
        else:
            # Values are both positive and negative, use symmetric normalization
            cmap = plt.get_cmap('coolwarm')
            norm = TwoSlopeNorm(vmin=-max_abs_value, vcenter=0, vmax=max_abs_value)
            normalized_values = values

        # Plot the lines
        fig, ax = plt.subplots()
        for (p1, p2), value in zip(
                [(p1, p2) for p1, p2, value in allsegments],
                normalized_values):
            x1, y1 = p1
            x2, y2 = p2
            ax.plot([x1, x2], [y1, y2], color=cmap(norm(value)), linewidth=2)

        # Set equal scaling
        ax.set_aspect('equal', adjustable='box')

        # Add a color bar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Value')

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Lines with Value-Based Coloring')
        plt.show()
        return
