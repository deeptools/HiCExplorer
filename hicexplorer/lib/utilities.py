import logging
log = logging.getLogger(__name__)


class Utilities():

    def __init__(self):
        pass

    def in_units(self, pBasePosition):
        pBasePosition = float(pBasePosition)
        log.debug("pBasePosition {}".format(pBasePosition))
        if pBasePosition > 1.5e6:
            labels = "{:.2f} ".format(pBasePosition / 1e6)
            labels += " Mbp"
        elif pBasePosition > 1500:
            labels = "{:.0f}".format(pBasePosition / 1e3)
            labels += " Kbp"
        else:
            labels = "{:.2f} ".format((pBasePosition))
            labels += " bp"
        return labels
    # def relabel_ticks(self, pXTicks):

    #     # log.debug('type pXTicks {} '.format(type(pXTicks)))
    #     if pXTicks[-1] > 1.5e6:
    #         labels = ["{:.2f} ".format(x / 1e6)
    #                 for x in pXTicks]
    #         labels[-2] += " Mbp"
    #     elif pXTicks[-1] > 1500:
    #         labels = ["{:.0f}".format(x / 1e3)
    #                 for x in pXTicks]
    #         labels[-2] += " Kbp"
    #     else:
    #         labels = ["{:.2f} ".format((x))
    #                 for x in pXTicks]
    #         labels[-2] += " bp"
    #     return labels
