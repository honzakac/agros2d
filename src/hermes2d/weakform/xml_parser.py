#!/usr/bin/python

from xml.dom import minidom
import os
from expression_parser import NumericStringParser
import module as md
from module import expression, weakform


class Config:     
    """Stores important data like paths, template names etc. """
    
    modules_dir = '../../../modules/'
    weakform_dir = './'
    factory_dir = '../' 
    doc_dir = './doc/'     
    templates = ['weakform_cpp.tem', 'weakform_h.tem', 'weakform_factory_h.tem']
    project_file = 'weakform.pri'


class XmlParser: 
    """ Manages modules xml files. """

    def __init__(self, modules):
        """ Reads xml files for required modules. """
        self.module_files = []
        for module in modules:
            self.module_files.append(module + '.xml')
                                       
        self.modules = []
        self.templates = dict()        
       
        # For generating c++ code are used templates containing structure of  
        # source and header files, which are in fact xml files.   
        for template_file in Config.templates:
            try:              
                fread = open(template_file, 'r')
                template = minidom.parse(fread).documentElement        
                fread.close()    
                self.templates[template_file] = template
            except IOError as (errno, strerror):
                print 'I/O error({0}): {1} "{2}".'.format(errno, strerror, 
														  template_file)
                raise
                                  
    def process(self): 
        """ Generates source and header c++ files. """                                          
        for module_file in self.module_files:        
            self.parse_xml_file(module_file)

        # removes pri file to makes qmake compiles new c++ files 
        try:
            os.remove(Config.weakform_dir + Config.project_file)
        except:
            pass
        
        files = []
        conditions = []

        for module in self.modules:                           
            module_files, module_conditions = module.get_code(self.templates)            
            module.write_code(Config.weakform_dir, self.templates)                                       
            conditions.extend(module_conditions)                                      
            files.extend(module_files)
        # writes weakform_factory.h         
        factory_code_str = ''    
        key = 'weakform_factory_h.tem'        
        node = self.templates[key].getElementsByTagName('head')[0]  
        factory_code_str += node.childNodes[0].nodeValue      
       
        for module_file in files:                    
            if module_file[::-1][:2][::-1] == '.h':                
                node = self.templates[key].getElementsByTagName('includes')[0]  
                string = node.childNodes[0].nodeValue 
                string = string.replace("general_weakform.h", module_file)
                factory_code_str += string

        weakform_temps = ['CustomMatrixFormVol','CustomVectorFormVol',
                           'CustomMatrixFormSurf','CustomVectorFormSurf', 'CustomEssentialFormSurf']                    
        for weakform_temp in weakform_temps:                            
            weakform_string = ''
            for condition in conditions:                            
                if condition[0] == weakform_temp:
                    weakform_string += condition[1]            
            node = self.templates[key].getElementsByTagName(weakform_temp)[0]                          
            string = node.childNodes[0].nodeValue      
            string = string.replace('//conditions', weakform_string)
            factory_code_str += string
                
        node = self.templates[key].getElementsByTagName('footer')[0]  
        factory_code_str += node.childNodes[0].nodeValue         
        factory_file = open(Config.factory_dir+'weakform_factory.h', 'w')
        factory_file.write(factory_code_str)
        factory_file.close()


    def gen_doc(self):
        index_string = '.. toctree::\n    :maxdepth: 2\n\n'
        for module in self.modules:
            index_string += '    ' + module.id + '.rst' + '\n'
            doc_file = open(Config.doc_dir + module.id + '.rst', 'w') 
            doc_file_str = module.gen_doc()
            doc_file.write(doc_file_str)
            doc_file.close()
        index_file = open(Config.doc_dir + 'index.rst', 'w')
        index_file.write(index_string)
        index_file.close()

    def parse_xml_file(self, filename):
        module_xml = md.parse(Config.modules_dir + filename)
        general = module_xml.get_general()
        module = Module()
        module.id = general.id
        module.name = general.name
        module.description = general.description

        constants_xml = module_xml.constants
        for const_xml in constants_xml.get_constant():
            const = Constant()
            const.id = const_xml.id
            const.value = float(const_xml.value)
            module.constants.append(const)
#
#        postprocessor_xml = module_xml.get_postprocessor()
#        variables_xml = postprocessor_xml.localvariables        
#        for variable_xml in variables_xml.get_localvariable():
#            variable = Variable()            
#            variable.id = variable_xml.id
#            variable.name = variable_xml.name
#            variable.short_name = variable_xml.shortname
#            variable.short_name_latex = variable_xml.shortname_latex
#            variable.short_name_html = variable_xml.shortname_html
#            variable.unit = variable_xml.unit
#            variable.unit_latex = variable_xml.unit_latex
#            variable.unit_html = variable_xml.unit_html
#            variable.type = variable_xml.type_
#            module.variables.append(variable)

        volume = Volume()
        volume_xml = module_xml.get_volume()        
        for quantity_xml in volume_xml.get_quantity():
            variable = Variable()
            variable.id = quantity_xml.id
            variable.short_name = quantity_xml.shortname
            module.variables.append(variable)
            variable.type = 'material'
        weak_forms_xml = volume_xml.get_weakforms()
        for weak_form_xml in weak_forms_xml.get_weakform():
            volume.name = weak_form_xml.analysistype
            for matrix_xml in weak_form_xml.matrix:
                for coor in ['planar', 'axi']:
                    weakform = WeakForm()
                    weakform.coordinate_type = coor
                    weakform.analysis_type = volume.name
                    if coor == 'planar':
                        weakform.expression = matrix_xml.planar
                    else:
                        weakform.expression = matrix_xml.axi
                    weakform.type = "matrix"
                    weakform.integral_type = 'vol'
                    weakform.i = matrix_xml.i
                    weakform.j = matrix_xml.j
                    volume.weakforms.append(weakform)
            for vector_xml in weak_form_xml.vector:
                for coor in ['planar', 'axi']:
                    weakform = WeakForm()
                    weakform.coordinate_type = coor
                    weakform.analysis_type = volume.name
                    if coor == 'planar':
                        weakform.expression = vector_xml.planar
                    else:
                        weakform.expression = vector_xml.axi
                    weakform.type = "vector"
                    weakform.integral_type = 'vol'
                    weakform.i = matrix_xml.i
                    weakform.j = matrix_xml.j
                    volume.weakforms.append(weakform)
            module.volumes.append(volume)            
            
            surface = Surface()
            surface_xml = module_xml.get_surface()
            quantities = []
            for quantity_xml in surface_xml.get_quantity():
                quantity = Variable()
                quantity.id = quantity_xml.id
                quantity.short_name = quantity_xml.shortname
                quantity.type = 'boundary'
                quantities.append(quantity)
            surface.variables = quantities
            weak_forms_xml = surface_xml.get_weakforms()
            for weakform_xml in weak_forms_xml.get_weakform():
                weakform = WeakForm()
                weakform.integral_type = 'surf'
                weakform.analysis_type = weakform_xml.analysistype
                weakform.default = weakform_xml.default
                for boundary_xml in weakform_xml.get_boundary():
                    boundary = Boundary()
                    boundary.weakforms = []
                    boundary.id = boundary_xml.id
                    boundary.name = boundary_xml.name
                    boundary.quantities = []
                    for quantity_xml in volume_xml.get_quantity():
                        quantity = Variable()
                        quantity.id = quantity_xml.id
                        quantity.short_name = quantity_xml.shortname
                        quantity.type = 'boundary'
                        boundary.quantities.append(quantity)
                    for matrix_xml in boundary_xml.matrix:
                        for coor in ['planar', 'axi']:
                            weakform = WeakForm()
                            weakform.coordinate_type = coor
                            weakform.analysis_type = volume.name
                            if coor == 'planar':
                                weakform.expression = matrix_xml.planar
                            else:
                                weakform.expression = matrix_xml.axi
                            weakform.type = "matrix"
                            weakform.integral_type = 'surf'
                            weakform.i = matrix_xml.i
                            weakform.j = matrix_xml.j
                            weakform.quantities = boundary.quantities
                            weakform.boundary_type = boundary_xml.id
                            surface.weakforms.append(weakform)
                        
                    for vector_xml in boundary_xml.vector:
                        for coor in ['planar', 'axi']:
                            weakform = WeakForm()
                            weakform.coordinate_type = coor
                            weakform.analysis_type = volume.name
                            if coor == 'planar':
                                weakform.expression = vector_xml.planar
                            else:
                                weakform.expression = vector_xml.axi
                            weakform.type = "vector"
                            weakform.integral_type = 'surf'
                            weakform.i = vector_xml.i
                            weakform.j = vector_xml.j
                            weakform.boundary_type = boundary_xml.id
                            surface.weakforms.append(weakform)
            module.surfaces.append(surface)
            self.modules.append(module)

class WeakForm:
    def __init__(self):        
        self.id = ''
        self.type = ''
        self.integral_type = ''        
        self.coordinate_type = ''
        self.analysis_type = ''
        self.boundary_type = ''        
        self.i = 0
        self.j = 0
        self.expression = ''
        self.name = ''
        self.variables = []
        
    def get_temp_class_name(self):
        class_name =  'Custom' + self.type.capitalize() + 'Form'  \
            + self.integral_type.capitalize()                        
        return class_name    
       
    def get_class_name(self):
        class_name =  'Custom' + self.type.capitalize() + 'Form'  \
            + self.integral_type.capitalize() + '_' + self.coordinate_type + '_' + str(self.i) \
            + '_'  + str(self.j)                        
        return class_name    
    
    def get_function_name(self):              
        function_name =  'custom_' + self.type + '_form_' + self.integral_type 
        return function_name    

    def get_factory_code(self, factory_template):
        if (self.type == 'essential'):
            string = factory_template.getElementsByTagName('condition_exact_solution')[0].childNodes[0].nodeValue                                
        elif (self.type == 'vector'):
            if(self.integral_type == 'vol'):
                string = factory_template.getElementsByTagName('condition_vector_vol')[0].childNodes[0].nodeValue
            else:    
                string = factory_template.getElementsByTagName('condition_vector_surf')[0].childNodes[0].nodeValue
        else:
            if(self.integral_type == 'vol'):
                string = factory_template.getElementsByTagName('condition_matrix_vol')[0].childNodes[0].nodeValue
            else:
                string = factory_template.getElementsByTagName('condition_matrix_surf')[0].childNodes[0].nodeValue
        
                     
        string = string.replace('class_name', self.id + '_'+  self.coordinate_type)
        string = string.replace('axi', 'axisymmetric')                        
        string = string.replace('row_index', str(self.i-1))                        
        string = string.replace('column_index', str(self.j-1))
        string = string.replace('boundary_type', self.boundary_type)        
        print string
        namespace = self.id.replace('_','')
        string = string.replace('namespace', namespace)
        function_name = self.get_class_name();        
        string = string.replace('FunctionName', function_name)                        
        factory_code = []
        factory_code.append(self.get_temp_class_name())
        factory_code.append(string)                        
        return factory_code
        
    def get_h_code(self, h_template):                                        
        h_code = ''   
        node = h_template.getElementsByTagName('variable_declaration')[0]
        variable_def_temp = node.childNodes[0].nodeValue                                                                                                                     
        for node in h_template.getElementsByTagName(self.get_function_name()):                        
            string = node.childNodes[0].nodeValue                                                                                                                                                                                  
            name = self.get_temp_class_name()             
            variable_defs = ''
            replaced_string = string.replace(name, name + '_' + self.coordinate_type     + '_' + str(self.i) \
                + '_'  + str(self.j))            
            for variable in self.variables:                                    
                variable_string = variable_def_temp.replace('variable_short', 
                                        variable.short_name)                    
                variable_string = variable_string.replace('variable', 
                                        variable.id)
                variable_defs += variable_string
                variable_defs = variable_defs.replace('volume', variable.type)                                                                                        
            replaced_string = replaced_string.replace('//variable_declaration', 
                                                              str(variable_defs))                             
            h_code += replaced_string                                                                                                                           
        return h_code 
    
    def get_cpp_code(self, cpp_template):                                
        function_types = ['','_value', '_ord', '_derivatives']
        cpp_code = ''        
        for function_type in function_types:         
            node = cpp_template.getElementsByTagName('variable_definition')[0]
            variable_def_temp = node.childNodes[0].nodeValue                                            
            for node in cpp_template.getElementsByTagName(self.get_function_name() + function_type):                        
                string = node.childNodes[0].nodeValue                                                                                                                                                                                                  
                name = self.get_temp_class_name()             
                replaced_string = string.replace(name, name + '_' + self.coordinate_type + '_' + str(self.i) \
                + '_'  + str(self.j))            
                if function_type == '':
                    variable_defs = '' ;                    
                    for variable in self.variables:                    
                        variable_string = variable_def_temp.replace('variable_short', 
                                        variable.short_name)                    
                        variable_string = variable_string.replace('variable', 
                                        variable.id)
                        variable_defs += variable_string
                        variable_defs = variable_defs.replace('material', variable.type)
                    replaced_string = replaced_string.replace('//variable_definition', 
                                                              str(variable_defs))     
                if function_type == '_ord':
                    expression = self.parse_expression(self.expression, True, '')                     
                else:
                    expression = self.parse_expression(self.expression, False, '')                                                                                                                                                                  
                if self.expression =='':               
                    replaced_string = ''
                else:                    
                    replaced_string = replaced_string.replace('//expression', 
                                expression) + '\n\n'                
                cpp_code += replaced_string                                                                                                                           
        return cpp_code            
            
    def parse_expression(self, expression, without_variables, output):
        replaces = { 'PI': 'M_PI',
                     'f': 'Util::scene()->problemInfo()->frequency',
                     'x': 'e->x[i]',
                     'y': 'e->y[i]',
                     'r': 'e->x[i]',
                     'z': 'e->y[i]',
                     'udr': 'u->dx[i]',
                     'vdr': 'v->dx[i]',
                     'udz': 'u->dy[i]',
                     'vdz': 'v->dy[i]',
                     'updr': 'u_ext[this->j]->dx[i]',
                     'updz': 'u_ext[this->j]->dy[i]',
                     'udx': 'u->dx[i]',
                     'vdx': 'v->dx[i]',
                     'udy': 'u->dy[i]',
                     'vdy': 'v->dy[i]',
                     'updx': 'u_ext[this->j]->dx[i]',
                     'updy': 'u_ext[this->j]->dy[i]',
                     'upval': 'u_ext[this->j]->val[i]',
                     'uval': 'u->val[i]',
                     'vval': 'v->val[i]', 
                     'uptval': 'ext->fn[this->i]->val[i]',
                     'deltat': 'Util::scene()->problemInfo()->timeStep.number()'
                     }
        
        
        latex_replaces = { '*': '\\cdot ',
                     'PI': '\\pi',
                     'EPS0': '\\varepsilon_0',
                     'epsr': '\\varepsilon_r',
                     'f': 'f',                     
                     'udx': '\\frac{\\partial u^{l}}{\\partial x}',
                     'udy': '\\frac{\\partial u^{l}}{\\partial y}',           
                     'udr': '\\frac{\\partial u^{l}}{\\partial r}',
                     'udz': '\\frac{\\partial u^{l}}{\\partial z}',
                     'vdx': '\\frac{\\partial v^{l}}{\\partial x}',
                     'vdy': '\\frac{\\partial v^{l}}{\\partial y}',           
                     'vdr': '\\frac{\\partial v^{l}}{\\partial r}',
                     'vdz': '\\frac{\\partial v^{l}}{\\partial z}',           
                     'uval': 'u',
                     'vval': 'v',
                     'upval': 'u^{l-1}',
                     'updx': '\\frac{\\partial u^{l-1}}{\\partial x}',
                     'updy': '\\frac{\\partial u^{l-1}}{\\partial y}',
                     'updr': '\\frac{\\partial u^{l-1}}{\\partial r}',
                     'updz': '\\frac{\\partial u^{l-1}}{\\partial z}',
                     'deltat': '\\delta t'                       
                     }            

                
        symbols = ['x', 'y', 'r', 'z', 'f', 'udr', 'udz', 'udx', 'udy',
                   'vdr', 'vdz', 'vdx', 'vdy', 'updr', 'updx', 'updy', 'updz',
                   'uval', 'vval', 'upval', 'deltat', 'uptval', 'PI']
                           
        variables = []
        variables_derivatives = []
        for variable in self.variables:
            symbols.append(variable.short_name)
            symbols.append("d" + variable.short_name)
            variables.append(variable.short_name)
            variables_derivatives.append("d" + variable.short_name)

        for const in self.constants:
            symbols.append(const.id)
        if output == 'latex':
            if not(expression.replace(' ','') == ''):
                parser = NumericStringParser(symbols, latex_replaces, variables,
                                             variables_derivatives,
                                             without_variables)
                expression_list = parser.parse(expression).asList()                                  
                parsed_exp = parser.translate_to_latex(expression_list)
            else:
                parsed_exp =''                             
        else:
            parser = NumericStringParser(symbols, replaces, variables, 
										 variables_derivatives,
										 without_variables)                        
            if not(expression.replace(' ','') == ''):
                expression_list = parser.parse(expression).asList()                                  
                parsed_exp = parser.translate_to_cpp(expression_list)                             
            else:
                parsed_exp =''
            parsed_exp = '(' + parsed_exp + ');'                                          
        return parsed_exp

class Volume:
    def __init__(self):
        self.id = ''        
        self.name = ''
        self.type = ''                
        self.weakforms = []        

class Surface:
    def __init__(self):
        self.id = ''
        self.name = ''
        self.type = ''        
        self.weakforms = []
        self.variables = []
                        
class Variable:
    def __init__(self):
        self.id = ''
        self.type = ''        
        self.name = ''
        self.short_name = ''
        self.units = ''  
        
    def write_cpp_code(self):
        pass

class Boundary:
    def __init__(self):
        pass

class Constant:
    def __init__(self):        
        self.id = ' '
        self.value = 0

class PartModule:
    def __init__(self):        
        self.id = ''
        self.name = ''
        self.description = ''        
        self.analysis = ''
        self.coordinate_type = '' 
        self.constants = []        
        self.weakforms = []        
        self.used_weakforms = set([])     
        self.forms_number = 0
        
class Module:
    def __init__(self):
        self.id = ''
        self.name = ''
        self.description = ''
        self.volumes = []
        self.surfaces = []
        self.constants = []
        self.variables = []

    def info(self):
        print 'ID: ', self.id
        print 'Name: ', self.name
        print 'Description: ', self.description
        i = 0
        print  '\nConstants:'
        print '--------------------------------'
        for constant in self.constants:
            print constant.id, constant.value
            i += 1                 
        print '\nAnalysis:'        
        print '--------------------------------'
        for volume in self.volumes:
            print '\n--------------------------------'
            print volume.name            
            print '--------------------------------'
            print '\nmatrix forms:'            
            for variable in self.variables:
                print variable.name
            for weakform in volume.weakforms:
                print (weakform.type, weakform.coordinate_type, 
                      weakform.integral_type)
                print weakform.i, weakform.j, weakform.expression

        print 'Boundaries'
        for boundary in self.boundaries:
            print 'Variables'            
            for quantity in boundary.quantities:
                print quantity[0], "   ", quantity[1]
            print 'Weakforms'        

    def underline(self, string, character):
        under_string = string + '\n' + character*len(string) + '\n'
        return under_string
    
    def gen_doc(self):
        doc_string = self.underline(self.name, '*')
        doc_string += self.underline('Description ', '=')\
                      + self.description + '\n\n'
        i = 0                                
        doc_string +=  self.underline('Constants', '=')        
        for constant in self.constants:
            doc_string += str(constant.id) + '  ' + str(constant.value) + '\n'
            i += 1                 
            doc_string += '\n' + self.underline('Analysis', '=') + '\n'        
        for volume in self.volumes:
            doc_string +='\n' + self.underline(volume.name.capitalize(),'-') + '\n'            
            doc_string += '\n' + self.underline('Domain weak forms:', '^') + '\n'            
#            for variable in self.variables:
#                doc_string += variable.name + '\n'
            
            for weakform in volume.weakforms:
                doc_string += weakform.type + weakform.coordinate_type 
                print weakform.coordinate_type
                doc_string += weakform.integral_type + '\n\n'
                doc_string += '.. math:: \n\n'                         
                doc_string += '    ' \
                              + weakform.parse_expression(weakform.expression, False,'latex') + '\n\n'                       
#        
        doc_string += '\n' + self.underline('Boundary conditions:', '-') + '\n'
        for boundary in self.boundaries:                                          
            doc_string += '\n' + self.underline('Variables:', '^') + '\n'
            for variable in boundary.variables:
                doc_string +=variable.short_name + '  ' + variable.name + '  [' + \
                              variable.unit + '] \n\n'                  
            doc_string += '\n' + self.underline('Weakforms', '^') + '\n'        
            for weakform in boundary.weakforms:
                doc_string += '.. math:: \n\n'                
                doc_string += '    ' + weakform.parse_expression(weakform.expression, False,'latex') + '\n\n'
        return doc_string
                                        
    def extract_modules(self):
        module_types = []
        part_modules = []
        part_module = PartModule()
        for volume in self.volumes:
            for weakform in volume.weakforms:
                part_module_id = self.id + '_' \
                    + weakform.analysis_type
                if (part_module_id in module_types):
                    index = module_types.index(part_module_id)
                    part_module = part_modules[index]
                else:
                    module_types.append(part_module_id)
                    part_module = PartModule()
                    part_module.name = self.name
                    part_module.id = part_module_id
                    part_module.description = self.description
                    part_module.coordinate_type = weakform.coordinate_type
                    part_module.constants = self.constants
                    part_module.volumes = self.volumes
                    part_module.analysis = volume.name
                    part_modules.append(part_module)
                weakform.variables = self.variables
                weakform.constants = self.constants
                weakform.id = self.id + '_' \
                    + weakform.analysis_type 
                part_module.weakforms.append(weakform)

        for surface in self.surfaces:
            for weakform in surface.weakforms:
                print weakform.id
                if weakform.id == 'vec':
                    pass
                part_module_id = self.id + '_' \
                    + weakform.analysis_type
                index = module_types.index(part_module_id)
                part_module = part_modules[index]
                weakform.id = self.id + \
                '_'+ weakform.analysis_type
                weakform.variables = surface.variables
                weakform.constants = self.constants
                part_module.weakforms.append(weakform)
#        for part_module in part_modules:
#             for weakform in part_module.weakforms:
#                 print weakform.id, weakform.i, weakform.j
#            for weakform in part_module.weakforms:
#                for variable in weakform.variables:
#                    print variable.short_name
        return part_modules;


    def get_code(self, param_templates):
        templates = dict() 
        templates['.cpp'] = param_templates['weakform_cpp.tem']
        templates['.h'] = param_templates['weakform_h.tem']
        file_strings = dict()
        part_modules = self.extract_modules()
        factory_codes = []
        for part_module in part_modules:
            filename = (part_module.id)
            for key in templates.iterkeys():
                file_string_name = filename + key
                node = templates[key].getElementsByTagName('head')[0]            
                string = node.childNodes[0].nodeValue                               
                file_strings[file_string_name] = string 
                node = templates[key].getElementsByTagName('includes')[0]            
                string = node.childNodes[0].nodeValue                                              
                string = string.replace('general_weakform', filename)           
                file_strings[file_string_name] += string  
                node = templates[key].getElementsByTagName('namespaces')[0]            
                string = node.childNodes[0].nodeValue                                              
                string = string.replace('general_weakform', filename)
                string = string.replace('_', '')
                file_strings[file_string_name] += string
                    
                class_names = set([])                            
               
                for weakform in part_module.weakforms:
                    if key == '.cpp':
                        class_names.add(weakform.get_class_name())
                        file_strings[file_string_name] += weakform.get_cpp_code(templates[key])            
                        factory_code =  weakform.get_factory_code(param_templates['weakform_factory_h.tem'])
                        factory_codes.append(factory_code)
                    if key == '.h':                        
                        file_strings[file_string_name] += weakform.get_h_code(templates[key])                                

                node = templates[key].getElementsByTagName('footer')[0]                            
                if key == '.cpp':                        
                    for class_name in class_names:
                        string = node.childNodes[0].nodeValue
                        string = string.replace('ClassName', class_name)                        
                        file_strings[file_string_name] += string             
               
                if key == '.h':       
                    string = node.childNodes[0].nodeValue                                     
                    file_strings[file_string_name] += string                                                                  

        return file_strings , factory_codes
       
                                       
    def write_code(self, weakform_dir, param_templates):                                                     
            weakform_pri_file = open(weakform_dir + 'weakform.pri', 'a')            
            files, conditions = self.get_code(param_templates)            

            for filename in files.iterkeys():                                            
                output_file = open(weakform_dir + filename , 'w')
                output_file.write(files[filename])            
                output_file.close()                            
               
                # append to weakform.pri                
                if filename[::-1][:4][::-1] == '.cpp':                
                    weakform_pri_file.write('SOURCES += ' + filename + '\n')            
#            weakform_pri_file.close()

def main():
    modules = ['acoustic']
    xml_parser = XmlParser(modules)
    xml_parser.process()
    #xml_parser.gen_doc()

if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    main()