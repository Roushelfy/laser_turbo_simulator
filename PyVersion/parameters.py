import math
import numpy as np
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GLU import gluPerspective
import cv2
from PIL import Image
from configparser import ConfigParser
config_parser = ConfigParser()
config_parser.read('parameters.cfg')
# config = config_parser['default']

# print(config['server'], config['user'], config['password'], 
    #   config['name'], config['charset'], sep='\n')

class Parameters():
    class ObjectFunction:
        def __init__(self):
            self.type = cfg['type']
            self.x_expr = cfg['x']
            self.y_expr = cfg['y']
            self.z_expr = cfg['z']

    def __init__(self):
        
        self.MAX_STOP_BREAK_r = float(config_parser['default']['max_time']) 
        self.CAPTURE_INTERVAL_r = float(config_parser['default']['capture_interval'])
        self.object_curve_function = self.object_curve_function(
            config_parser['Object_Function']['type'],
            config_parser['Object_Function']['x_expr'],
            config_parser['Object_Function']['y_expr'],
            config_parser['Object_Function']['z_expr']
        )
        self.object_update_frequency = float(config_parser['velocity']['object_update_frequency'])
        self.laser_update_frequency = float(config_parser['velocity']['laser_update_frequency'])
        self.sensors = self.generate_sensors(
            np.array((0.0, 0.0, 0.0), dtype=float),
            config_parser['sensors']['type'],
            int(config_parser['sensors']['number']),
            float(config_parser['sensors']['radius'])
        )
        self.laser = self.Laser_property(
            float(config_parser['laser']['wavelength']),
            float(config_parser['laser']['waist_radius']),
            float(config_parser['laser']['scale_factor'])
        )
        self.laser_move_function = self.laser_move_function(float(config_parser['laser']['move_alpha']))
        # init_laser_orientation = 0.0, 0.0, 1.0 in config_parser
        self.init_laser_position = np.array(eval(config_parser['laser']['init_position']), dtype=float)
        self.init_laser_orientation = np.array(eval(config_parser['laser']['init_orientation']), dtype=float)
        self.object_max_radius = 0.05

    class Laser_property():
        def __init__(self, wavelength:float, waist_radius:float, scale_factor:float):
            self.wavelength = wavelength
            self.waist_radius = waist_radius
            self.scale_factor = scale_factor #intensity 


    def object_curve_function(self, type:str, x_expr:str, y_expr:str, z_expr:str):
        def Lfunc(t:float):
            x = 0 + 1 * math.sin(t)
            y = 0 + 0.5 * math.sin(t)
            z = 6.00
            p_object = np.array((x, y, z), dtype=float)
            return np.array((x, y, z), dtype=float)

        def Cfunc(t:float):
            x = -0.3 + 0.3 * math.cos(t)
            y = 0.0 + 0.8 * math.sin(t)
            z = 6.00
            p_object = np.array((x, y, z), dtype=float)
            print(x,y,z)
            return np.array((x, y, z), dtype=float)

        def userDefine(t:float):
            namespace = {'math': math, 't': t}
            x = eval(x_expr, namespace)
            y = eval(y_expr, namespace)
            z = eval(z_expr, namespace)
            return np.array((x, y, z), dtype=float)
        print('type = ',type)
        print('x_expr = ',x_expr)
        print('y_expr = ',y_expr)
        print('z_expr = ',z_expr)

        if type == 'L':
            return Lfunc
        elif type == 'C':
            return Cfunc
        elif type == 'userDefine':
            return userDefine
        else:
            raise Exception('Unknown object function type')

    def generate_sensors(self,center:np.array, type:str, number:int, radius:float):

        def circle_generator(center:np.array, number:int, radius:float):
            for i in range(number):
                x = radius * math.cos(i * 2 * math.pi / number)
                y = radius * math.sin(i * 2 * math.pi / number)
                z = 0.0
                yield np.array((x, y, z), dtype=float) + center

        if type == 'Circle':
            sensors = []
            for sensor in circle_generator(center, number, radius):
                sensors += [sensor]
        return np.array(sensors)

    def laser_move_function(self,alpha):
        def move_func1(intensity,sensors,v_orientation):
            rel_v = (max(intensity) - min(intensity)) / min(intensity)
            # rate = float(math.exp(rel_v) ** 2) # 差距越大调整幅度越大
            rate = 1
            max3_index = np.argsort(intensity)[-3:]
            vec_sum = np.array((0.0,0.0,0.0),dtype=float)
            for i in max3_index:
                vec_sum += sensors[i]
            vec_sum /= 3
            # print('offset = ', self.object_position - self.laser_orientation / self.laser_orientation[2] * 6.0)
            print('vec_sum = ',vec_sum * alpha )
            return v_orientation + alpha * rate * vec_sum
        return move_func1


class VideoSaver():
    
    def render_ray(self,p_object, start1 ,end1, start2, end2):
        def draw_ray(start, end, color):
            glBegin(GL_LINES)
            glColor3fv(color)
            glVertex3fv(start)
            glVertex3fv(end)
            glEnd()

        def draw_point(position, color, size):
            # print('draw',position)
            glEnable(GL_POINT_SMOOTH)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

            glColor3fv(color)
            glPointSize(size)

            glBegin(GL_POINTS)
            glVertex3fv(position)
            glEnd()

            glDisable(GL_POINT_SMOOTH)
            glDisable(GL_BLEND)

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        left, right, bottom, top, near_clip, far_clip = -1, 1, -1, 1, 0.5, 20.0
        glOrtho(left, right, bottom, top, near_clip, far_clip)
        glTranslatef(0.0, 0.0, -5)
        
        try:
            draw_ray(start1, end1, (0.5, 0.5, 0.5))  # Blue color for the first ray
            draw_ray(start2, end2, (0.68, 0.85, 0.1))  # Light blue color for the second ray
        except:
            pass
        draw_point(p_object, (1, 0, 0), 10.0)  # Red point for the object

        # output x,y,z 3 axis
        draw_ray((0, 0, 0), (10, 0, 0), (0, 0, 1))
        draw_ray((0, 0, 0), (0, 10, 0), (0, 0, 1))
        draw_ray((0, 0, 0), (0, 0, 10), (0, 0, 1))
        glRotatef(45, 0, 0, 1)


    def __init__(self):
        pygame.init()
        display = (800, 600)
        pygame.display.set_mode(display, DOUBLEBUF | OPENGL)
        glEnable(GL_DEPTH_TEST)
        fourcc = cv2.VideoWriter_fourcc(*'XVID')
        self.out = cv2.VideoWriter('output.avi', fourcc, 20.0, (800, 600))

    def append_frame(self, laser_position, laser_orientation, object_position):
        start1 = laser_position
        distence = np.linalg.norm(object_position - laser_position)
        end1 = laser_position + laser_orientation * distence
        start2 = (2 * object_position - laser_position) * 0.5
        end2 = start2 - laser_orientation * distence
        self.render_ray(object_position, start1, end1, start2, end2)
        glReadBuffer(GL_FRONT) 
        data = glReadPixels(0, 0, 800, 600, GL_RGB, GL_UNSIGNED_BYTE)
        image = Image.frombytes("RGB", (800, 600), data)
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        #save to test.png
        # image.save('test.png', 'PNG')
        self.out.write(np.array(image))
        
class Simulation():
    def __init__(self):
        self.parameters = Parameters()
        self.move_period = max(1,self.parameters.laser_update_frequency / self.parameters.object_update_frequency)
        self.MAX_STOP_BREAK_i = int(self.parameters.MAX_STOP_BREAK_r / self.parameters.object_update_frequency)
        self.CAPTURE_INTERVAL_i = int(self.parameters.CAPTURE_INTERVAL_r / self.parameters.object_update_frequency)
        # print('move_period = ',self.move_period)
        self.object_position = self.parameters.object_curve_function(0)
        self.laser_position = self.parameters.init_laser_position
        self.laser_orientation = self.parameters.init_laser_orientation
        self.video_saver = VideoSaver()

    def move_object(self,t):
        self.object_position = self.parameters.object_curve_function(t)
    
    def query_intensity(self):
        def distance_to_line(point, pstart, orientation):
            assert np.linalg.norm(orientation) != 0
            return np.linalg.norm(np.cross(point - pstart, orientation)) / np.linalg.norm(orientation)

        def gaussian_beam_intensity(dest, start, orientation, wavelength, waist_radius):
            distance = np.linalg.norm(dest - start)
            r = distance_to_line(dest, start, orientation)
            w = waist_radius * math.sqrt(1 + (distance / (math.pi * waist_radius ** 2 / wavelength)) ** 2)
            intensity = math.exp(-2 * r ** 2 / w ** 2) * self.parameters.laser.scale_factor
            return intensity

        r = distance_to_line(self.object_position, self.laser_position, self.laser_orientation)
        if r > self.parameters.object_max_radius:
            return [0 for i in range(len(self.parameters.sensors))]

        new_laser_start = (2 * self.object_position - self.laser_position) * 0.5
        new_laser_orientation = self.laser_orientation * -1
        new_laser_orientation = new_laser_orientation / np.linalg.norm(new_laser_orientation)

        results = []
        for sensor in self.parameters.sensors:
            strength = gaussian_beam_intensity(sensor, new_laser_start, new_laser_orientation,
                self.parameters.laser.wavelength, self.parameters.laser.waist_radius)
            results += [strength]
        return np.array(results)


    def move_laser(self):
        intensity = self.query_intensity()
        self.laser_orientation = self.parameters.laser_move_function(intensity, self.parameters.sensors, self.laser_orientation)
        self.laser_position = self.laser_orientation

    def start(self):
        t = 0
        while True:
            self.move_object(t * self.parameters.object_update_frequency)
            if t % self.move_period == 0:
                self.move_laser()
            if t % self.CAPTURE_INTERVAL_i == 0:
                self.video_saver.append_frame(self.laser_position, self.laser_orientation, self.object_position)
            t += 1
            # print(self.laser_orientation)
            if t == self.MAX_STOP_BREAK_i:
                break


if __name__ == "__main__":
    simulation = Simulation()
    simulation.start()