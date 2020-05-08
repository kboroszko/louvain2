//
// Created by kajetan on 08.05.2020.
//

#include "gui.h"

#include <gtk/gtk.h>
#include <pthread.h>

int value;

void *argument_thread(void *args){
    while(value < 1000){
        value++;
        sleep(1);
    }
}

static void print_hello(GtkWidget *widget, gpointer data) {
    g_print("Hello World\n");
    GtkTextBuffer *buffer = (GtkTextBuffer *) data;
    char str[100];
    sprintf(str, "value:%d\n", value);
    gtk_text_buffer_set_text (buffer, str, -1);
}

static void activate(GtkApplication *app, gpointer user_data) {
    GtkWidget *window;
    GtkWidget *grid;
    GtkWidget *button;
//    GtkWidget *button_box;
    GtkWidget *text_box;
    GtkTextBuffer **buffer = (GtkTextBuffer **) user_data;

    window = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW (window), "Window");
    gtk_window_set_default_size(GTK_WINDOW (window), 400, 600);

    /* Here we construct the container that is going pack our buttons */
    grid = gtk_grid_new ();
    /* Pack the container in the window */
    gtk_container_add (GTK_CONTAINER (window), grid);


//    button_box = gtk_button_box_new(GTK_ORIENTATION_HORIZONTAL);
//    gtk_container_add(GTK_CONTAINER (window), button_box);

    text_box = gtk_text_view_new();
    *buffer = gtk_text_view_get_buffer (GTK_TEXT_VIEW (text_box));
    gtk_text_buffer_set_text (*buffer, "Hello, this is some text", -1);

//    gtk_container_add(GTK_CONTAINER (window), button_box);
//    gtk_container_add(GTK_CONTAINER (window), text_box);

    button = gtk_button_new_with_label("Hello World");
    g_signal_connect (button, "clicked", G_CALLBACK(print_hello), (void*)*buffer);
//    g_signal_connect_swapped (button, "clicked", G_CALLBACK(gtk_widget_destroy), window);
//    gtk_container_add(GTK_CONTAINER (button_box), button);
    gtk_grid_attach (GTK_GRID (grid), button, 0, 0, 1, 1);

    gtk_grid_attach (GTK_GRID (grid), text_box, 0, 1, 2, 10);

    gtk_widget_show_all(window);
}

int main(int argc, char **argv) {
    GtkApplication *app;
    GtkTextBuffer * text_buf;
    int status;
    pthread_t hello_world_thread;

    int result = pthread_create(&hello_world_thread, NULL, argument_thread, NULL);
    if (result != 0) {
        perror("Could not create thread.");
    }

    app = gtk_application_new("org.gtk.example", G_APPLICATION_FLAGS_NONE);
    g_signal_connect (app, "activate", G_CALLBACK(activate), &text_buf);


    status = g_application_run(G_APPLICATION (app), argc, argv);



    g_object_unref(app);

    return status;
}