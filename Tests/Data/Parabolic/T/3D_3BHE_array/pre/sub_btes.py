import numpy as np
import pandas as pd

from tespy import subsys, cmp, con
from tespy.helpers import MyComponentError


class btes_para(subsys.subsystem):

    def __init__(self, label, num_btes, **kwargs):

        if not isinstance(label, str):
            raise MyComponentError('Subsystem label must be of type str!')

        elif len([x for x in [';', ', ', '.'] if x in label]) > 0:
            raise MyComponentError('Can\'t use ' + str([';', ', ', '.']) + ' ',
                                   'in label.')
        else:
            self.label = label

        if num_btes <= 1:
            raise MyComponentError('Minimum number of compressors is 2.')
        else:
            self.num_btes = num_btes

        # set default values
        for key in self.attr():
            self.__dict__.update({key: np.nan})
            self.__dict__.update({key + '_set': False})

        self.subsys_init()
        self.set_attr(**kwargs)
        

    def attr(self):

        values = ([n for n in subsys.subsystem.attr(self) if
                   n != 'num_i' and n != 'num_o'])

        for i in range(self.num_btes):
            j = str(i)
            # BTES
            values += ['Q' + j, 'pr' + j, 'T_out' + j, 'T_in' + j ,'L_bhe'+ j,
                       'ks_bhe'+ j,'D_bhe'+ j,'hydro_group']

#            # pipe feed flow (from BTES)
#            values += ['pr_pf' + j, 'Q_pf' + j, 'dT_pf' + j,
#                       'L_pf' + j, 'ks_pf' + j, 'D_pf' + j]
#            
#            # pipe back flow (to BTES)
#            values += ['pr_pb' + j, 'Q_pb' + j, 'dT_pb' + j,
#                       'L_pb' + j, 'ks_pb' + j, 'D_pb' + j]
            
#        # ambient temperature
#        values += ['t_a', 't_a_design']
#        
#        # subsurface temperature
#        values +=[ 't_sub', 't_sub_design']
            
        return values

    def create_comps(self):

        self.num_i = 1
        self.num_o = 1
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)

#        self.vessel = []
        self.btes = []
#        self.pipe_feed = []
#        self.pipe_back = []
        
        self.splitter = cmp.splitter(label=self.label + '_splitter',
                                     num_out=self.num_btes)
        self.merge = cmp.merge(label=self.label + '_merge',
                               num_in=self.num_btes)

        for i in range(self.num_btes):
            j = str(i +1)
            self.btes += [cmp.heat_exchanger_simple(label=self.label +
                                                    '_' + j,
                                                    offdesign=['zeta', 'kA'])]
#            self.pipe_feed += [cmp.pipe(label=self.label +
#                                        '_pipe feed_' + j,
#                                        offdesign=['zeta', 'kA'])]
#            self.pipe_back += [cmp.pipe(label=self.label +
#                                        '_pipe back_' + j,
#                                        offdesign=['zeta', 'kA'])]

    def set_comps(self):

        for i in range(self.num_btes):
            j = str(i)
            self.btes[i].set_attr(pr=self.get_attr('pr' + j),
                                  Q=self.get_attr('Q' + j),
                                  D=self.get_attr('D_bhe' + j),
                                  L=self.get_attr('L_bhe' + j),
                                  ks=self.get_attr('ks_bhe' + j),
                                  hydro_group=self.get_attr('hydro_group'))
#                                  t_a_design=self.get_attr('t_sub_design'),
#                                  t_a=self.get_attr('t_sub'))
            
#            self.pipe_feed[i].set_attr(pr=self.get_attr('pr_pf' + j),
#                                       Q=self.get_attr('Q_pf' + j),
#                                       D=self.get_attr('D_pf' + j),
#                                       L=self.get_attr('L_pf' + j),
#                                       ks=self.get_attr('ks_pf' + j))
##                                       t_a_design=self.get_attr('t_a_design'),
##                                       t_a=self.get_attr('t_a'))
#            
#            self.pipe_back[i].set_attr(pr=self.get_attr('pr_pb' + j),
#                                       Q=self.get_attr('Q_pb' + j),
#                                       D=self.get_attr('D_pb' + j),
#                                       L=self.get_attr('L_pb' + j),
#                                       ks=self.get_attr('ks_pb' + j))
##                                       t_a_design=self.get_attr('t_a_design'),
##                                       t_a=self.get_attr('t_a'))

    def create_conns(self):

        self.conns = []
        
        self.conns += [con.connection(self.inlet, 'out1',
                                      self.splitter, 'in1')]
        self.conns += [con.connection(self.merge, 'out1',
                                      self.outlet, 'in1')]

        for i in range(self.num_btes):
            j = str(i + 1)
#            if i > 0:
#                # mass flow distributes equally to all pipes
#                self.conns += [con.connection(self.splitter, 'out' + j,
#                                              self.pipe_back[i], 'in1',
#                                              m=con.ref(self.conns[0],
#                                                        1 / self.num_btes, 0))]
#            else:
#                self.conns += [con.connection(self.splitter, 'out' + j,
#                                              self.pipe_back[i], 'in1')]
#                
#            self.conns += [con.connection(self.splitter, 'out' + j,
#                                          self.pipe_back[i], 'in1')]
#                
#            self.conns += [con.connection(self.pipe_back[i], 'out1',
#                                          self.btes[i], 'in1',
#                                          design=['T'])]
#            self.conns += [con.connection(self.btes[i], 'out1',
#                                          self.pipe_feed[i], 'in1',
#                                          design=['T'])]
#            self.conns += [con.connection(self.pipe_feed[i], 'out1',
#                                          self.merge, 'in' + j)]
                                         
            self.conns += [con.connection(self.splitter, 'out' + j,
                                          self.btes[i], 'in1')]
            self.conns += [con.connection(self.btes[i], 'out1',
                                          self.merge, 'in' + j)]
    def set_conns(self):

        if not hasattr(self, 'nw'):
            self.create_network()

        i = 0
        for bt in self.btes:
            j = str(i)
   
            inc = pd.DataFrame()
            inc['t'] = self.nw.conns.t == bt
            inc['t_id'] = self.nw.conns.t_id == 'in1'
            c, cid = inc['t'] == True, inc['t_id'] == True
            inc = inc.index[c & cid][0]

            outc = pd.DataFrame()
            outc['s'] = self.nw.conns.s == bt
            outc['s_id'] = self.nw.conns.s_id == 'out1'
            c, cid = outc['s'] == True, outc['s_id'] == True
            outc = outc.index[c & cid][0]

            inc.set_attr(T=self.get_attr('T_in' + j))
            outc.set_attr(T=self.get_attr('T_out' + j))
            i += 1
            
#        i = 0
#        for pipe in self.pipe_feed:
#            j = str(i)
#
#            inc = pd.DataFrame()
#            inc['t'] = self.nw.conns.t == pipe
#            inc['t_id'] = self.nw.conns.t_id == 'in1'
#            c, cid = inc['t'] == True, inc['t_id'] == True
#            inc = inc.index[c & cid][0]
#
#            outc = pd.DataFrame()
#            outc['s'] = self.nw.conns.s == pipe
#            outc['s_id'] = self.nw.conns.s_id == 'out1'
#            c, cid = outc['s'] == True, outc['s_id'] == True
#            outc = outc.index[c & cid][0]
#
#            if self.get_attr('dT_pf' + j + '_set'):
#                inc.set_attr(T=con.ref(outc, 1,
#                                       self.get_attr('dT_pf' + j)))
#            else:
#                inc.set_attr(T=np.nan)
#            i += 1
#
#        i = 0
#        for pipe in self.pipe_back:
#            j = str(i)
#
#            inc = pd.DataFrame()
#            inc['t'] = self.nw.conns.t == pipe
#            inc['t_id'] = self.nw.conns.t_id == 'in1'
#            c, cid = inc['t'] == True, inc['t_id'] == True
#            inc = inc.index[c & cid][0]
#
#            outc = pd.DataFrame()
#            outc['s'] = self.nw.conns.s == pipe
#            outc['s_id'] = self.nw.conns.s_id == 'out1'
#            c, cid = outc['s'] == True, outc['s_id'] == True
#            outc = outc.index[c & cid][0]
#
#            if self.get_attr('dT_pb' + j + '_set'):
#                outc.set_attr(T=con.ref(inc, 1,
#                              -self.get_attr('dT_pb' + j)))
#            else:
#                outc.set_attr(T=np.nan)
#            i += 1
#            