using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using OpenTK.Graphics.OpenGL;

namespace WindowsFormsApp9
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            //glControl1.MouseWheel += new System.Windows.Forms.MouseEventHandler(glControl1_MouseWhell);
            glControl1.MouseWheel += glControl1_MouseWhell;
            glControl1.KeyPress += glControl1_KeyPress;

        }
        double k = 1;
        double q = 1.2;
        int l = 0;
        double deltaX = 0;
        double deltaY = 0;
        private void glControl1_MouseWhell(object sender, MouseEventArgs e)
        {

            //if (e.Delta > 0)
            //{


            //        trackBar6.Value++;

            //}
            //else
            //{

            //    trackBar6.Value--; 

            //}
            l += e.Delta / 120;

            label3.Text = l.ToString();


        }
        private void glControl1_KeyPress(object sender, KeyPressEventArgs e)
        {

            if (e.KeyChar==(char)Keys.S)
            {
         
                deltaY -= 0.1/k;
                e.Handled = true;

            }
            if (e.KeyChar == (char)Keys.W)
            {

                deltaY += 0.1/k;
                e.Handled = true;

            }
            if (e.KeyChar == (char)Keys.A)
            {
                deltaX-= 0.1/k;
                e.Handled = true;
            }
            if (e.KeyChar == (char)Keys.D)
            {
                deltaX += 0.1/k;
                e.Handled = true;
            }
           
        }


        private void Form1_Load(object sender, EventArgs e)
        {
            GL.ClearColor(1.0f, 1.0f, 1.0f, 1.0F);
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();

            GL.Ortho(-5, 5, -5, 5, 5, -5);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();

        }


        void Draw()
        {
            GL.Clear(ClearBufferMask.ColorBufferBit);
            drowall();
            GL.PushMatrix();
           
            GL.Scale(k, k, 1);
            GL.Translate(-deltaX, -deltaY, 1);
            //  Dif();

            // DrawCurve();


            DrawPhas();
            PointDraw();
             PhasmasVectors();

            GL.PopMatrix();

            glControl1.SwapBuffers();
        }
        private void drowall()
        {
            timer1.Enabled = true;
            double hx = 0.1;
            double hy = 0.1;

            GL.PointSize(10.0f);
            GL.LineWidth(2.0f);

            GL.Begin(PrimitiveType.Lines);
            GL.Color3(0.8f, 0.8f, 0.8f);
            for (double l = -5; l <= 5; l += 0.5) //шаг 10 часть от 2
            {
                GL.Vertex2(-5, l);
                GL.Vertex2(5, l);
            }
            for (double l = -5; l <= 5; l += 0.5)
            {
                GL.Vertex2(l, -5);
                GL.Vertex2(l, 5);
            }




            GL.End();
            GL.Begin(PrimitiveType.Lines);
            GL.Color3(0.4f, 0.4f, 0.4f);
            for (double l = -5; l <= 5; l += 1) //шаг половина 2
            {
                GL.Vertex2(-5, l);
                GL.Vertex2(5, l);
            }
            for (double l = -5; l <= 5; l += 1)
            {
                GL.Vertex2(l, -5);
                GL.Vertex2(l, 5);
            }



            GL.End();
            GL.Begin(PrimitiveType.Lines);
            GL.Color3(0.2f, 0.2f, 0.2f);
            for (double l = -5; l <= 5; l += 2) //шаг 2
            {
                GL.Vertex2(-5, l);
                GL.Vertex2(5, l);
            }
            for (double l = -5; l <= 5; l += 2)
            {
                GL.Vertex2(l, -5);
                GL.Vertex2(l, 5);
            }



            GL.End();

            GL.Begin(PrimitiveType.Lines);
            GL.Color3(0.0f, 0.0f, 0.0f);
            GL.Vertex2(5, 0);
            GL.Vertex2(-5, 0);
            GL.Vertex2(0, 5);
            GL.Vertex2(0, -5);
            GL.End();

            GL.Begin(PrimitiveType.Triangles);

            GL.Color3(0.0f, 0.0f, 0.0f);
            GL.Vertex2(5, 0);
            GL.Vertex2(5 - 3 * hx, 3 * hy);
            GL.Vertex2(5 - 3 * hx, -3 * hy);
            GL.Vertex2(0, 5);
            GL.Vertex2(3 * hx, 5 - 3 * hy);
            GL.Vertex2(-3 * hx, 5 - 3 * hy);
            GL.End();
            double h1 = 2.0;
            double h2 = 1.0;
            double h3 = 0.5;

            GL.Begin(PrimitiveType.Lines);
            for (double l = -5.0; l <= 5 - hx; l += h3)
            {
                GL.Vertex2(l, 0.15);
                GL.Vertex2(l, -0.15);
            }
            for (double l = -5.0; l <= 5 - hx; l += h2)
            {
                GL.Vertex2(l, 0.15);
                GL.Vertex2(l, -0.15);
            }
            for (double l = -5.0; l <= 5 - hx; l += h3)
            {
                GL.Vertex2(0.25, l);
                GL.Vertex2(-0.25, l);
            }
            for (double l = -5.0; l <= 5 - hx; l += h2)
            {
                GL.Vertex2(0.25, l);
                GL.Vertex2(-0.25, l);
            }
            for (double l = -5.0; l <= 5 - hx; l += h1)
            {
                GL.Vertex2(0.5, l);
                GL.Vertex2(-0.5, l);
            }
            for (double l = -5.0; l <= 5 - hx; l += h1)
            {
                GL.Vertex2(0.5, l);
                GL.Vertex2(-0.5, l);
            }
            for (double l = -5.0; l <= 5 - hx; l += 0.1)
            {
                GL.Vertex2(0.05, l);
                GL.Vertex2(-0.05, l);
            }
            for (double l = -5.0; l <= 5 - hx; l += 0.1)
            {
                GL.Vertex2(l, 0.05);
                GL.Vertex2(l, -0.05);
            }
            GL.End();


        }

        int ValuePresent = 0;
        int ValuePresent2 = 0;
        int Value3 = 0;
        int Value4 = 0;
        void PointDraw()
        {
            GL.PointSize(4.0f);
            GL.Color3(Color.Blue);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex2(deltaX, deltaY);
            GL.End();

        }


        //private void AprSmas()
        //{
        //    AprS = new double[N + 1, 2];

        //    AprS[0, 0] = ValuePresent/20.0;
        //    AprS[0, 1] = ValuePresent2/20.0;



        //}
        private void trackBar1_Scroll(object sender, EventArgs e)
        {




        }
        int N = 1;
        private void timer1_Tick(object sender, EventArgs e)
        {

            ValuePresent = trackBar1.Value;
            ValuePresent2 = trackBar2.Value;
            Value3 = trackBar4.Value;
            Value4 = trackBar5.Value;

            N = trackBar3.Value;
            textBox1.Text = ValuePresent.ToString();
            textBox3.Text = N.ToString();
            textBox2.Text = ValuePresent2.ToString();
            textBox4.Text = Value3.ToString();
            textBox5.Text = Value4.ToString();

            k = Math.Pow(q, l);
            label7.Text = (1/k* -5).ToString();
            label8.Text = (1/k * 5).ToString();
            label9.Text = (1/k * -5).ToString();
            label10.Text = (1/k * 5).ToString();
            textBox6.Text = (deltaX/k).ToString();
            textBox7.Text = (deltaY/k). ToString();
     
            Phas();
          
            Draw();

        }
        /// <summary>
        ///масташаб
        ///вектор дельта 
        ///расчет
        ///транслайт
        ///в 0
        ///обратно масштибирование 
        ///.в началае scale потом translate 
        ///
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void button1_Click(object sender, EventArgs e)
        {
            Draw();
        }

        private void label1_Click(object sender, EventArgs e)
        {

        }
        double f1(double x, double y)
        {
            return Value3/100.0 *y+ValuePresent/100.0 * x ;
        }
        double f2(double x, double y)
        {
            return y * ValuePresent2/100.0 + Value4/100.0*x;
        }

        //double[,] AprS ;
        //double X1(double t,double C1,double C2)
        //{
        //    return (3*C1*Math.Sin(6*t)+3*C2*Math.Cos(6*t)) * Math.Exp(-2*t);
        //}
        //double X2(double t, double C1, double C2)
        //{
        //    return ((2*C2- C1) * Math.Sin(6 * t) -(C2+2*C1) * Math.Cos(6 * t)) * Math.Exp(-2 * t);
        //}

        //void Dif()
        //{


        //    double dt = 10.0/ N; //T=10;N=tracbar;
        //    for (int i = 1; i < N +1; i++)
        //    {
        //        AprS[i, 0] = AprS[i - 1, 0] + dt * f1(AprS[i - 1, 0], AprS[i - 1, 1]);
        //        AprS[i, 1] = AprS[i - 1, 1] + dt * f2(AprS[i - 1, 0], AprS[i - 1, 1]);


        //    }
        //    GL.LineWidth(1.0f);
        //    GL.Color3(Color.Red);
        //    GL.Begin(PrimitiveType.LineStrip);
        //    for(int i = 0; i < N + 1; i++)
        //    {
        //        GL.Vertex2(AprS[i, 0], AprS[i, 1]);
        //    }
        //    GL.End();
        //}

    
        double[,] PhasmasF(double a, double b)
        {
            double dt = 10.0 / N;
            double[,] Phasmas = new double[N + 1, 2];
            Phasmas[0, 0] = a;
            Phasmas[0, 1] = b;
            for (int i = 1; i < N + 1; i++)

            {

                Phasmas[i, 0] = Phasmas[i - 1, 0] +dt * f1(Phasmas[i - 1, 0], Phasmas[i - 1, 1]);
                Phasmas[i, 1] = Phasmas[i - 1, 1] + dt * f2(Phasmas[i - 1, 0], Phasmas[i - 1, 1]);



            }
            return Phasmas;

        }

        double[,] PhasmasFV(double a, double b)
        {
            double dt = 10.0 / 1;
            double[,] Phasmas = new double[10 + 1, 2];
            Phasmas[0, 0] = a;
            Phasmas[0, 1] = b;
            for (int i = 1; i < 10 + 1; i++)

            {

                Phasmas[i, 0] = Phasmas[i - 1, 0] + dt * f1(Phasmas[i - 1, 0], Phasmas[i - 1, 1]);
                Phasmas[i, 1] = Phasmas[i - 1, 1] + dt * f2(Phasmas[i - 1, 0], Phasmas[i - 1, 1]);



            }

            return Phasmas;

        }
        List<double[,]> Phase;
      

    
        void DrawPhas()
        {
            GL.LineWidth(1.0f);
            GL.Color3(Color.Red);

            for (int i = 0; i < Phase.Count; i++)
            {
                GL.Begin(PrimitiveType.LineStrip);
                for (int j = 0; j < N + 1; j++)
                {
                    
                    GL.Vertex2(Phase[i][j, 0], Phase[i][j, 1]);
                    

                }
                GL.End();
            }

        }
        void Phas()
        {
            Phase = new List<double[,]>();
            for (double i = -5; i < 5; i += 0.5)

            {
                for (double j = -5; j < 5; j += 0.5)
                {
                    Phase.Add(PhasmasF(i/ k+deltaX, j / k+deltaY));

                }
            }


        }


      
    

        void DrowVector(double a, double b, double h, float L, double min, double max)
        {

            float colormax = 1.0f;
            float colormin = 1.0f;

            double H1 = h * f1(a, b) / Math.Sqrt(f1(a, b) * f1(a, b) + f2(a, b) * f2(a, b));
            double H2 = h * f2(a, b) / Math.Sqrt(f1(a, b) * f1(a, b) + f2(a, b) * f2(a, b));
            GL.LineWidth(2.0f);



            GL.Begin(PrimitiveType.Lines);

            GL.Color3(L * colormax + (1 - L) * 0, L * 0 + colormin * (1 - L), L * 0.0f + 0.0f * (1 - L));

            GL.Vertex2(a, b);
            GL.Vertex2(a - H1, b - H2);
            GL.Vertex2(a, b);
            GL.Vertex2(a - H1 + 0.1, b - H2 + 0.1);
            GL.Vertex2(a, b);
            GL.Vertex2(a - H1 - 0.1, b - H2 + 0.1);
            GL.Vertex2(a, b - 0.1);
            GL.End();


        }

        void PhasmasVectors()
        {
            List<double[,]> Phase = new List<double[,]>();

            for (double i = -5; i < 6; i += 1)

            {
                for (double j = -5; j < 6; j += 1)
                {
                    Phase.Add(PhasmasFV(i, j));

                }
            }
            double[] color = new double[Phase.Count - 1];

            double max = 0;

            double min = Math.Sqrt(f1(Phase[0][0, 0], Phase[0][0, 1]) * f1(Phase[0][0, 0], Phase[0][0, 1]) + f2(Phase[0][1, 0], Phase[0][1, 1]) * f2(Phase[0][1, 0], Phase[0][1, 1])); ;
            for (int i = 0; i < Phase.Count - 1; i++)
            {

                {

                    if (Math.Sqrt(f1(Phase[i][0, 0], Phase[i][0, 1]) * f1(Phase[i][0, 0], Phase[i][0, 1]) + f2(Phase[i][1, 0], Phase[i][1, 1]) * f2(Phase[i][1, 0], Phase[i][1, 1])) > max)
                    {
                        max = Math.Sqrt(f1(Phase[i][0, 0], Phase[i][0, 1]) * f1(Phase[i][0, 0], Phase[i][0, 1]) + f2(Phase[i][1, 0], Phase[i][1, 1]) * f2(Phase[i][1, 0], Phase[i][1, 1]));
                    }

                    if (Math.Sqrt(f1(Phase[i][1, 0], Phase[i][1, 1]) * f1(Phase[i][1, 0], Phase[i][1, 1]) + f2(Phase[i][1, 0], Phase[i][1, 1]) * f2(Phase[i][1, 0], Phase[i][1, 1])) < min)
                    {
                        min = Math.Sqrt(f1(Phase[i][0, 0], Phase[i][0, 1]) * f1(Phase[i][0, 0], Phase[i][0, 1]) + f2(Phase[i][1, 0], Phase[i][1, 1]) * f2(Phase[i][1, 0], Phase[i][1, 1]));
                    }
                }
            }
            float lam = 0;
            for (int i = 0; i < Phase.Count - 1; i++)

            {

                lam = (float)((Math.Sqrt(f1(Phase[i][0, 0], Phase[i][0, 1]) * f1(Phase[i][0, 0], Phase[i][0, 1]) + f2(Phase[i][1, 0], Phase[i][1, 1]) * f2(Phase[i][1, 0], Phase[i][1, 1])) - min) / (max - min));




                DrowVector(Phase[i][0, 0], Phase[i][0, 1], 0.4, lam, min, max);



            }
        }
        private void trackBar2_Scroll(object sender, EventArgs e)
        {

        }

        private void trackBar3_Scroll(object sender, EventArgs e)
        {

        }

        private void glControl1_Load(object sender, EventArgs e)
        {

        }

        private void trackBar6_Scroll(object sender, EventArgs e)
        {

        }



        private void textBox6_TextChanged(object sender, EventArgs e)
        {

        }

        

        private void glControl1_MouseDown(object sender, MouseEventArgs e)
        {

        }

        private void glControl1_KeyDown(object sender, KeyEventArgs e)
        {

        }

        private void label6_Click(object sender, EventArgs e)
        {

        }
    }
}
